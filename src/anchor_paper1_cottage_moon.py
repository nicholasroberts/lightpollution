"""Anchor accepted Cottage RAW camera maps to north using the clicked observed Moon."""
import argparse,json,shutil,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from reconstruct_spider_eye_contrast import digest

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--review-root',type=Path,required=True)
    p.add_argument('--points',type=Path,default=ROOT/'config/paper1/Cottage_reference_moon_points.json')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    points=json.loads(a.points.read_text());source=a.review_root/'registration'
    lp=source/'locked_site_references.json'
    if points['site']!='Cottage' or digest(lp)!=points['source_lock_sha256']:
        p.error('Moon selections refer to a different source lock')
    doc=json.loads(lp.read_text());lock=doc['sites']['Cottage']
    if lock['status']!='numerically_locked_awaiting_visual_review':p.error('Registration pending')
    record=next(r for r in lock['raw_inputs'] if Path(r['path']).stem=='1')
    if digest(record['path'])!=points['reference_sha256'] or record['sha256']!=points['reference_sha256']:
        p.error('Reference image changed')
    raw=Path(points['raw_reference_path'])
    if digest(raw)!=record['source_raw_sha256']:p.error('Original NEF changed')
    import exifread,ephem
    with raw.open('rb') as stream:tags=exifread.process_file(stream,details=False)
    meta=points['raw_metadata']
    for key in ['EXIF DateTimeOriginal','EXIF OffsetTimeOriginal']:
        if str(tags.get(key,''))!=meta[key]:p.error('Original RAW timestamp/offset changed')
    local=datetime.strptime(meta['EXIF DateTimeOriginal'],'%Y:%m:%d %H:%M:%S')
    utc=datetime.fromisoformat(local.isoformat()+meta['EXIF OffsetTimeOriginal']).astimezone(timezone.utc)
    if utc.isoformat()!='2018-08-21T20:23:49+00:00':p.error('Unexpected Cottage capture UTC')
    selected=[r for r in points['points'] if r['label'].casefold()=='moon']
    if len(selected)!=1:p.error('Exactly one point labelled Moon is required')
    moon_xy=np.asarray(selected[0]['pixel_xy'],float);center=np.asarray(lock['reference_center_xy_px'])
    radius=lock['rim_radius_px']
    implied=90-np.rad2deg(2*np.arcsin(np.linalg.norm(moon_xy-center)/(np.sqrt(2)*radius)))
    observer=ephem.Observer();observer.lat='55.426041';observer.lon='14.231644'
    observer.elevation=0;observer.pressure=0;observer.date=utc.strftime('%Y/%m/%d %H:%M:%S')
    moon=ephem.Moon(observer);sun=ephem.Sun(observer)
    elevation=float(np.rad2deg(float(moon.alt)));azimuth=float(np.rad2deg(float(moon.az)))
    if elevation<=0 or abs(implied-elevation)>5:p.error('Moon elevation/radius inconsistent with the reviewed geometry')
    anchor=dict(site_id=lock['site_id'],local_time=meta['EXIF DateTimeOriginal'],
        camera_utc_offset=meta['EXIF OffsetTimeOriginal'],timezone='EXIF explicit offset; not inferred from Sweden',
        utc_offset_hours=datetime.fromisoformat(local.isoformat()+meta['EXIF OffsetTimeOriginal']).utcoffset().total_seconds()/3600,
        utc=utc.isoformat().replace('+00:00','Z'),timestamp_source='Original NEF DateTimeOriginal and OffsetTimeOriginal',
        reference_field_date='2018-08-21',moon_elevation_deg=elevation,moon_compass_azimuth_deg=azimuth,
        moon_illuminated_percent=float(moon.phase),sun_elevation_deg=float(np.rad2deg(float(sun.alt))),
        sun_compass_azimuth_deg=float(np.rad2deg(float(sun.az))),
        ephemeris='PyEphem '+ephem.__version__+'; pressure=0; height=0; reference exposure start, whole seconds')
    registration=dict(method='Observed Moon bearing and fixed ephemeris azimuth; inherited zenith-axis equisolid projection',
        observed_moon_pixel_xy=moon_xy.tolist(),image_implied_moon_elevation_deg=float(implied),
        elevation_residual_deg=float(implied-elevation),moon_radius_adjusted=False,
        handedness='reversed',parity_support='P1 near Vega and P4 near Arcturus under Moon-anchored orientation; no parity fit to polarization',
        star_points=points['points'][1:],source_points_sha256=digest(a.points),
        source_lock_sha256=digest(lp),polarization_used_for_orientation=False,
        limitations=['Moon bearing alone does not calibrate camera tilt or parity.',
            'Moon radial/elevation discrepancy retained, not fitted away.',
            'Reference timestamp is exposure start, not the effective midpoint of the 12-minute sequence.'])
    lock.update(astronomical_anchor=anchor,moon_pixel_xy=moon_xy.tolist(),
        moon_position_status='observed_reference_Moon_clicked_by_investigator',
        celestial_orientation=dict(handedness='reversed',moon_registration=registration),
        north_orientation_status='observed_Moon_anchored_RAW_timestamp_explicit_offset')
    a.output.mkdir(parents=True)
    dest=a.output/'registration';shutil.copytree(source,dest)
    outlock=dest/'locked_site_references.json';outlock.write_text(json.dumps(doc,indent=2)+'\n')
    coords=dest/'Cottage/reference_coords.json'
    if coords.exists():coords.write_text(json.dumps(lock,indent=2)+'\n')
    cp=dest/'reference_sha256.json';checks=json.loads(cp.read_text())
    checks[outlock.name]=digest(outlock)
    if coords.exists():checks[str(coords.relative_to(dest))]=digest(coords)
    cp.write_text(json.dumps(checks,indent=2)+'\n')
    (a.output/'moon_north_anchor.json').write_text(json.dumps(dict(astronomical_anchor=anchor,registration=registration),indent=2)+'\n')
    src=Path(__file__).resolve().parent
    def run(name,*args):subprocess.run([sys.executable,str(src/name),*map(str,args)],check=True)
    # Re-audit against the derived lock; no registration or exposure normalization.
    run('audit_spider_measured_stokes.py','--reference-lock',outlock,'--sites','Cottage',
        '--raw-manifest-root',Path(record['path']).parent,'--output',a.output/'audit')
    run('rotate_spider_measured_north.py','--reference-lock',outlock,
        '--audit',a.output/'audit/Cottage','--site','Cottage','--handedness','reversed',
        '--output',a.output/'north')
    run('plot_paper1_measured_panel.py','--reference-lock',outlock,'--audit',a.output/'audit/Cottage',
        '--north',a.output/'north','--site','Cottage','--output',a.output/'panel')
    print(json.dumps(anchor,indent=2));print('Moon radial discrepancy (degrees):',implied-elevation)


if __name__=='__main__':main()
