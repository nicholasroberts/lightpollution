"""Prepare derived PennCC stellar orientation and confirmed WCpitch lock for Figure 1D."""
import argparse, hashlib, json, shutil
from datetime import datetime, timezone
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def stellar_solution(points,observer,center,radius):
    import ephem
    xy=np.array([p['pixel_xy'] for p in points['points']],float)
    delta=xy-center;beta=np.arctan2(delta[:,1],delta[:,0])
    image_el=np.pi/2-2*np.arcsin(np.linalg.norm(delta,axis=1)/(radius*np.sqrt(2)))
    stars=[]
    for point in points['points']:
        name=points['star_identities'][point['label']]
        star=ephem.star(name);star.compute(observer)
        stars.append(dict(label=point['label'],name=name,pixel_xy=point['pixel_xy'],
            elevation_deg=float(np.rad2deg(float(star.alt))),azimuth_deg=float(np.rad2deg(float(star.az))),
            magnitude=float(star.mag)))
    alt=np.deg2rad([s['elevation_deg'] for s in stars])
    az=np.deg2rad([s['azimuth_deg'] for s in stars])
    solutions=[]
    for sign in [-1,1]:
        offset=float(np.angle(np.mean(np.exp(1j*(az-sign*beta)))))
        predicted=sign*beta+offset
        dots=np.sin(image_el)*np.sin(alt)+np.cos(image_el)*np.cos(alt)*np.cos(predicted-az)
        errors=np.rad2deg(np.arccos(np.clip(dots,-1,1)))
        solutions.append(dict(azimuth_vs_pixel_angle_sign=sign,
            azimuth_offset_deg=float(np.rad2deg(offset)),
            rms_position_deg=float(np.sqrt(np.mean(errors**2))),
            point_residuals_deg=errors.tolist(),
            handedness='reversed' if sign==-1 else 'original'))
    solution=min(solutions,key=lambda s:s['rms_position_deg'])
    if solution['rms_position_deg']>3:raise ValueError('Stellar residual exceeds reviewed 3-degree check')
    for star,el in zip(stars,image_el):star['image_projection_elevation_deg']=float(np.rad2deg(el))
    return solution,solutions,stars


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--registration-root',type=Path,required=True)
    p.add_argument('--metadata-audit',type=Path,required=True)
    p.add_argument('--points',type=Path,default=ROOT/'config/paper1/PennCC_reference_star_points.json')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    import ephem
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import tifffile
    points=json.loads(a.points.read_text())
    source=a.registration_root/'runs/PennCC'
    lp=source/'locked_site_references.json'
    if digest(lp)!=points['source_lock_sha256']:p.error('Point selections refer to a different source lock')
    penn=json.loads(lp.read_text())['sites']['PennCC']
    ref=next(r for r in penn['raw_inputs'] if Path(r['path']).stem=='1')
    if ref['sha256']!=points['reference_sha256'] or digest(ref['path'])!=points['reference_sha256']:
        p.error('Reference TIFF differs from picked image')
    center=np.asarray(penn['reference_center_xy_px'],float);radius=penn['rim_radius_px']
    if not np.allclose(center,points['reference_center_xy_px']) or radius!=points['rim_radius_px']:
        p.error('Reference projection differs')
    rows=[r for r in json.loads(a.metadata_audit.read_text()) if Path(r['path']).parts[-2]=='PennCC']
    meta=min(rows,key=lambda r:r['metadata']['EXIF DateTimeOriginal'])['metadata']
    local=datetime.strptime(meta['EXIF DateTimeOriginal'],'%Y:%m:%d %H:%M:%S')
    if local.strftime('%Y-%m-%d %H:%M:%S')!='2019-04-20 22:24:20':
        p.error('PennCC timestamp differs from investigator-confirmed capture')
    utc=datetime.fromisoformat(local.isoformat()+meta['EXIF OffsetTimeOriginal']).astimezone(timezone.utc)
    observer=ephem.Observer();observer.lat='52.558887';observer.lon='-2.145724'
    observer.elevation=0;observer.pressure=0;observer.date=utc.strftime('%Y/%m/%d %H:%M:%S')
    solution,alternatives,stars=stellar_solution(points,observer,center,radius)
    moon=ephem.Moon(observer);alt=float(moon.alt);az=float(moon.az)
    offset=np.deg2rad(solution['azimuth_offset_deg']);sign=solution['azimuth_vs_pixel_angle_sign']
    beta=(az-offset)/sign
    xy=center+radius*np.sqrt(1-np.sin(alt))*np.array([np.cos(beta),np.sin(beta)])
    anchor=dict(site_id='s08',local_time=meta['EXIF DateTimeOriginal'],timezone='Europe/London',
        utc_offset_hours=datetime.fromisoformat(local.isoformat()+meta['EXIF OffsetTimeOriginal']).utcoffset().total_seconds()/3600,
        utc=utc.isoformat().replace('+00:00','Z'),timestamp_source='NEF DateTimeOriginal and explicit OffsetTimeOriginal',
        reference_field_date='2019-04-20',moon_elevation_deg=float(np.rad2deg(alt)),
        moon_compass_azimuth_deg=float(np.rad2deg(az)),moon_illuminated_percent=float(moon.phase),
        flags=['Moon pixel predicted from accepted four-star orientation; Moon is not directly visible'],
        ephemeris='PyEphem '+ephem.__version__+'; geometric topocentric; pressure=0; height=0')
    provenance=dict(method='Four identified stars; fixed centre and equisolid radius; unweighted circular-mean bearing offset; compare both parities',
        solution=solution,parity_alternatives=alternatives,stars=stars,
        source_points_sha256=digest(a.points),source_lock_sha256=digest(lp),
        source_metadata_audit_sha256=digest(a.metadata_audit),predicted_moon_pixel_xy=xy.tolist(),
        polarization_used=False,centre_radius_or_tilt_fitted=False,investigator_accepted='2026-11-24',
        limitation='Residual includes projection/centroid/level-camera error; no independent holdout star')
    wcsource=a.registration_root/'runs/WCpitch'
    wc=json.loads((wcsource/'locked_site_references.json').read_text())['sites']['WCpitch']
    if not np.allclose(wc.get('moon_pixel_xy'),[621,638]):p.error('WCpitch candidate differs from investigator-confirmed point')
    a.output.mkdir(parents=True)
    for site,original in [('PennCC',source),('WCpitch',wcsource)]:
        dest=a.output/'runs'/site;shutil.copytree(original,dest)
        lockpath=dest/'locked_site_references.json'
        allsites=json.loads(lockpath.read_text());lock=allsites['sites'][site]
        if site=='PennCC':
            lock.update(astronomical_anchor=anchor,moon_pixel_xy=xy.tolist(),
                moon_position_status='predicted_from_reviewed_stars_and_ephemeris_not_observed',
                celestial_orientation=dict(handedness=solution['handedness'],stellar_registration=provenance),
                north_orientation_status='investigator_accepted_four_star_orientation')
        else:
            lock['moon_position_status']='observed_candidate_confirmed_by_investigator_2026-11-24'
        lockpath.write_text(json.dumps(allsites,indent=2)+'\n')
        cp=dest/'reference_sha256.json';checks=json.loads(cp.read_text())
        checks[lockpath.name]=digest(lockpath);cp.write_text(json.dumps(checks,indent=2)+'\n')
    report=dict(site='PennCC',astronomical_anchor=anchor,**provenance)
    (a.output/'PennCC_stellar_registration.json').write_text(json.dumps(report,indent=2)+'\n')
    image=tifffile.imread(ref['path']);lo,hi=np.nanpercentile(image,[5,99.9])
    fig,ax=plt.subplots(figsize=(10,7));ax.imshow(image,cmap='gray',vmin=lo,vmax=hi)
    for star in stars:
        aa=np.deg2rad(star['azimuth_deg']);ee=np.deg2rad(star['elevation_deg'])
        bb=(aa-offset)/sign
        predicted=center+radius*np.sqrt(1-np.sin(ee))*np.array([np.cos(bb),np.sin(bb)])
        ax.plot(*star['pixel_xy'],'rx');ax.plot(*predicted,'c+',ms=9)
        ax.annotate(star['label']+' '+star['name'],star['pixel_xy'],xytext=(10,10),textcoords='offset points',color='red')
    ax.plot(*xy,'y+',ms=12);ax.annotate('Predicted Moon (not observed)',xy,xytext=(-240,-30),textcoords='offset points',color='yellow')
    ax.set_title('PennCC: red observed stars; cyan predictions; yellow predicted Moon\nStellar RMS '+format(solution['rms_position_deg'],'.2f')+' degrees; no polarization fit')
    fig.savefig(a.output/'PennCC_stellar_registration.png',dpi=200);plt.close(fig)
    print(json.dumps(report,indent=2));print('Derived registration root:',a.output)


if __name__=='__main__':main()
