"""Accepted ChurchHill star orientation and zoology Moon: fits and Supp Figure 3."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace
import numpy as np
from correct_spider_analyzer_sign import apply as sign_correct, digest
from prepare_paper1_stellar_registration import stellar_solution
from run_paper1_oxford_natural_fits import export_reused_stokes, render

SITES=('ChurchHill','zoology')
EXPECTED={
'ChurchHill':'b5e38fae962ca026a26786c698633a035c1e01a7077fc0449d59f9acf806b754',
'zoology':'c81554432b3fc143f270afab9ab6b130aba58c90cb9ca53e8dd375834a7d3839'}

def run(a):
    import ephem
    if a.output.exists() or a.figure_copy.exists():
        raise ValueError('Use fresh output and figure-copy directories')
    points=json.loads(a.points.read_text())
    if points['source_lock_sha256']!=EXPECTED['ChurchHill']:
        raise ValueError('Star picks do not match reviewed ChurchHill lock')
    for site in SITES:
        path=a.registration_root/'runs'/site/'locked_site_references.json'
        if digest(path)!=EXPECTED[site]:raise ValueError('Reviewed source lock changed: '+site)
        lock=json.loads(path.read_text())['sites'][site]
        if lock.get('radiometry_status')!='black_subtracted_linear_camera_green':
            raise ValueError('Require linear camera green: '+site)
    source=json.loads((a.registration_root/'runs/ChurchHill/locked_site_references.json').read_text())['sites']['ChurchHill']
    ref=next(r for r in source['raw_inputs'] if Path(r['path']).stem=='1')
    if digest(ref['path'])!=points['reference_sha256'] or ref['sha256']!=points['reference_sha256']:
        raise ValueError('Picked ChurchHill reference changed')
    if source['astronomical_anchor']['utc']!=points['utc']:
        raise ValueError('ChurchHill capture UTC differs')
    c=np.array(source['reference_center_xy_px']);radius=source['rim_radius_px']
    if not np.allclose(c,points['reference_center_xy_px']) or radius!=points['rim_radius_px']:
        raise ValueError('Star projection changed')
    observer=ephem.Observer();observer.lat=str(points['latitude_deg']);observer.lon=str(points['longitude_deg'])
    observer.date=points['utc'].replace('T',' ');observer.pressure=0;observer.elevation=0
    solution,alternatives,stars=stellar_solution(points,observer,c,radius)
    if solution['handedness']!='reversed' or abs(solution['rms_position_deg']-points['solution']['rms_position_deg'])>.05:
        raise ValueError('Recomputed stellar solution differs from accepted result')
    moon=ephem.Moon(observer);alt=float(moon.alt);az=float(moon.az)
    beta=(az-np.deg2rad(solution['azimuth_offset_deg']))/solution['azimuth_vs_pixel_angle_sign']
    moonxy=c+radius*np.sqrt(1-np.sin(alt))*np.array([np.cos(beta),np.sin(beta)])
    a.output.mkdir(parents=True)
    reports={}
    def call(name,*args):
        subprocess.run([sys.executable,str(Path(__file__).with_name(name)),*map(str,args)],check=True)
    for site in SITES:
        base=a.output/site;original=a.registration_root/'runs'/site/'locked_site_references.json'
        # Existing helper independently determines each site's nominal analyzer sign.
        sign_correct(SimpleNamespace(reference_lock=original,sites=[site],output=base/'registration'))
        lockpath=base/'registration/locked_site_references.json'
        payload=json.loads(lockpath.read_text());lock=payload['sites'][site]
        expected_sign=-1 if site=='ChurchHill' else 1
        if lock['analyzer_stokes_convention']['applied_U_multiplier']!=expected_sign:
            raise ValueError('Unexpected analyzer sign: '+site)
        lock['visual_alignment_status']='investigator_accepted_2026-11-27'
        if site=='ChurchHill':
            lock['moon_pixel_xy']=moonxy.tolist()
            lock['moon_position_status']='predicted_from_reviewed_stars_and_ephemeris_not_observed'
            lock['north_orientation_status']='investigator_accepted_four_star_orientation'
            lock['celestial_orientation']=dict(handedness=solution['handedness'],
                stellar_registration=dict(solution=solution,alternatives=alternatives,stars=stars,
                    points_sha256=digest(a.points),investigator_accepted='2026-11-27',
                    polarization_used=False))
            lock['astronomical_anchor'].update(moon_elevation_deg=float(np.rad2deg(alt)),
                moon_compass_azimuth_deg=float(np.rad2deg(az)),moon_illuminated_percent=float(moon.phase),
                flags=['Moon predicted from accepted four-star orientation; not observed'],
                ephemeris='PyEphem '+ephem.__version__+'; pressure=0; height=0')
        else:
            if not np.allclose(lock['moon_pixel_xy'],[1420,631]):
                raise ValueError('Zoology Moon differs from accepted point')
            lock['moon_position_status']='observed_candidate_confirmed_by_investigator_2026-11-27'
            lock['celestial_orientation']=dict(handedness='reversed',
                note='Retain campaign camera-to-sky parity convention; not selected by polarization fit')
        lockpath.write_text(json.dumps(payload,indent=2)+'\n')
        (base/'registration'/site/'reference_coords.json').write_text(json.dumps(lock,indent=2)+'\n')
        checks={str(p.relative_to(lockpath.parent)):digest(p) for p in lockpath.parent.rglob('*')
                if p.is_file() and p.name!='reference_sha256.json'}
        (lockpath.parent/'reference_sha256.json').write_text(json.dumps(checks,indent=2)+'\n')
        stokes=lockpath.parent/site/'registration/linear_stokes.npz'
        with np.load(stokes,allow_pickle=False) as z:arrays={k:np.array(z[k]) for k in ('I','Q','U')}
        export_reused_stokes(base/'camera_export',site,stokes,lockpath,arrays)
        (base/'source_lock_path.txt').write_text(str(lockpath.resolve())+'\n')
        call('rotate_spider_measured_north.py','--reference-lock',lockpath,
             '--audit',base/'camera_export','--site',site,'--handedness','reversed','--output',base/'north')
        call('fit_paper1_solar_lunar_angle.py','--reference-lock',lockpath,
             '--north',base/'north','--site',site,'--neutral-separation-deg',70,
             '--fit-neutral-separation','--output',base/'fit')
        reports[site]=json.loads((base/'fit/solar_lunar_fit.json').read_text())
    summary=[]
    for site,fit in reports.items():
        row=dict(site=site,utc=fit['solar']['utc'],
            sun_elevation_deg=fit['solar']['sun_elevation_deg'],
            sun_compass_azimuth_deg=fit['solar']['sun_compass_azimuth_deg'],
            solar_fraction=fit['solar_fraction'],lunar_fraction=fit['lunar_fraction'],
            neutral_separation_deg=fit['neutral_separation_deg'],
            fraction_on_boundary=fit['fraction_on_boundary'],separation_on_boundary=fit['separation_on_boundary'])
        for model in ('lunar_only','solar_only','combined'):
            row[model+'_withheld_rms_deg']=fit['metrics'][model]['pooled_withheld']['rms_deg']
        summary.append(row)
    (a.output/'site_fit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (a.output/'site_fit_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    render(a.output,a.figure_copy,sites=SITES,
           title='Supplementary Figure 3 — ChurchHill and zoology\nShared mixture/separation per site; no ALAN or fitted dilution',
           stem='supp_fig03')
    for site in SITES:
        dest=a.figure_copy/'results'/site
        shutil.copytree(a.output/site/'fit',dest)
    for name in ('site_fit_summary.json','site_fit_summary.csv'):
        shutil.copy2(a.output/name,a.figure_copy/name)
    manifest={str(p.relative_to(a.figure_copy)):digest(p) for p in a.figure_copy.rglob('*') if p.is_file()}
    (a.figure_copy/'archive_sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(summary,indent=2));print('Saved',a.figure_copy)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--registration-root',type=Path,default=Path('data/processed/spider_registration_green_all_sites_01'))
    p.add_argument('--points',type=Path,default=Path('config/paper1/ChurchHill_reference_star_points.json'))
    p.add_argument('--output',type=Path,default=Path('data/processed/paper1_natural_fit_supp03_01'))
    p.add_argument('--figure-copy',type=Path,default=Path('figures/paper1/supplemental/supp_fig03/v01'))
    a=p.parse_args()
    try:run(a)
    except (OSError,ValueError,KeyError) as exc:p.error(str(exc))

if __name__=='__main__':main()
