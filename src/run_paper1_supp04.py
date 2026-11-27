"""Fit the confirmed stannes/opensquare Moon-anchored patterns; Supp Figure 4."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace
import numpy as np
from correct_spider_analyzer_sign import apply as sign_correct, digest, evidence
from run_paper1_oxford_natural_fits import export_reused_stokes, render

SITES=('stannes','opensquare')

def preflight(root,review):
    for site in SITES:
        path=root/'runs'/site/'locked_site_references.json'
        accepted=review[site]
        if digest(path)!=accepted['lock_sha256']:
            raise ValueError('Reviewed source lock changed: '+site)
        lock=json.loads(path.read_text())['sites'][site]
        if lock.get('radiometry_status')!='black_subtracted_linear_camera_green':
            raise ValueError('Require linear green: '+site)
        if not np.array_equal(lock['moon_pixel_xy'],accepted['moon_pixel_xy']):
            raise ValueError('Confirmed Moon differs: '+site)
        if lock['astronomical_anchor']['utc']!=accepted['utc']:
            raise ValueError('Reviewed capture UTC differs: '+site)
        if evidence(lock)[0]!=accepted['analyzer_sequence_sign']:
            raise ValueError('Signed analyzer evidence differs: '+site)
        checks=json.loads((path.parent/'reference_sha256.json').read_text())
        for rel,sha in checks.items():
            source=path.parent/rel
            if path.parent.resolve() not in source.resolve().parents or digest(source)!=sha:
                raise ValueError('Source checksum differs: '+str(source))
        stokes=path.parent/site/'registration/linear_stokes.npz'
        if digest(stokes)!=accepted['stokes_sha256']:
            raise ValueError('Reviewed Stokes changed: '+site)

def run(a):
    if a.output.exists() or a.figure_copy.exists():
        raise ValueError('Use fresh output and figure-copy directories')
    review=json.loads(a.review.read_text())
    preflight(a.registration_root,review)
    reports={}
    def call(name,*args):
        subprocess.run([sys.executable,str(Path(__file__).with_name(name)),*map(str,args)],check=True)
    for site in SITES:
        base=a.output/site
        original=a.registration_root/'runs'/site/'locked_site_references.json'
        sign_correct(SimpleNamespace(reference_lock=original,sites=[site],output=base/'registration'))
        lockpath=base/'registration/locked_site_references.json'
        payload=json.loads(lockpath.read_text());lock=payload['sites'][site]
        lock['moon_position_status']='observed_candidate_confirmed_by_investigator_2026-11-27'
        # Moon acceptance does not imply a fresh review of every exposure overlay.
        lock['moon_localization']['status']='investigator_confirmed_image_coordinate'
        lock['celestial_orientation']=dict(handedness='reversed',
            note='Retained campaign camera-to-sky parity convention; not selected by fit')
        anchor=lock['astronomical_anchor']
        anchor['flags']=[f for f in anchor.get('flags',[]) if not f.startswith('Ephemeris-only forward anchor:')]
        anchor['flags'].append('Observed Moon coordinate confirmed by investigator 2026-11-27')
        lockpath.write_text(json.dumps(payload,indent=2)+'\n')
        (lockpath.parent/site/'reference_coords.json').write_text(json.dumps(lock,indent=2)+'\n')
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
           title="Supplementary Figure 4 — St Anne's and Open Square\nShared mixture/separation per site; no ALAN or fitted dilution",
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
    p.add_argument('--review',type=Path,default=Path('config/paper1/supp04_reviewed_moons.json'))
    p.add_argument('--output',type=Path,default=Path('data/processed/paper1_natural_fit_supp04_01'))
    p.add_argument('--figure-copy',type=Path,default=Path('figures/paper1/supplemental/supp_fig03/v01'))
    a=p.parse_args()
    try:run(a)
    except (OSError,ValueError,KeyError) as exc:p.error(str(exc))

if __name__=='__main__':main()
