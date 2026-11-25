"""Run existing audited Stokes/north/angle-fit stages for Figure 1D; no registration refit."""
import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from paper1_figure_archive import archive_figure, digest

SITES={'PennCC':'2019-04-20','WCpitch':'2019-05-14'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--registration-root',type=Path)
    p.add_argument('--raw-manifest-root',type=Path)
    p.add_argument('--metadata-audit',type=Path)
    p.add_argument('--reuse-root',type=Path,help='Previously completed 1D run with verified SITE/north and source_lock_path.txt')
    p.add_argument('--figure-copy',type=Path,required=True)
    p.add_argument('--replace-figure-copy',action='store_true')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    if a.figure_copy.exists() and not a.replace_figure_copy:p.error('Figure copy exists; use explicit --replace-figure-copy')
    if a.reuse_root is None and not any([a.registration_root,a.raw_manifest_root,a.metadata_audit]):
        candidates={}
        for path in Path('data/processed').glob('*/PennCC/source_lock_path.txt'):
            root=path.parents[1]
            try:
                signature=[]
                for site in SITES:
                    lock=Path((root/site/'source_lock_path.txt').read_text().strip())
                    nr=json.loads((root/site/'north/north_rotation.json').read_text())
                    if nr['site']!=site or nr['source_lock_sha256']!=digest(lock):raise ValueError('Lineage mismatch')
                    signature.extend([nr['source_lock_sha256'],nr['source_arrays_sha256'],digest(root/site/'north/north_stokes.npz')])
                candidates.setdefault(tuple(signature),root)
            except (OSError,ValueError,KeyError):continue
        if len(candidates)!=1:p.error('Cannot uniquely select previous verified 1D products; supply --reuse-root PATH_TO_PREVIOUS_1D_RUN')
        a.reuse_root=next(iter(candidates.values()))
        print('Reusing verified north products:',a.reuse_root,flush=True)
    if a.reuse_root is not None:
        if any([a.registration_root,a.raw_manifest_root,a.metadata_audit]):p.error('Use either --reuse-root or original reconstruction inputs')
        a.output.mkdir(parents=True)
        src=Path(__file__).resolve().parent
        def run(script,*args):subprocess.run([sys.executable,str(src/script),*map(str,args)],check=True)
        for site in SITES:
            old=a.reuse_root/site;lock=Path((old/'source_lock_path.txt').read_text().strip())
            nr=json.loads((old/'north/north_rotation.json').read_text())
            if nr['site']!=site or nr['source_lock_sha256']!=digest(lock):p.error('Reused north/lock mismatch: '+site)
            base=a.output/site;base.mkdir()
            (base/'source_lock_path.txt').write_text(str(lock.resolve())+'\n')
            # Reuse the accepted celestial/Stokes product without another resampling.
            (base/'north').symlink_to((old/'north').resolve(),target_is_directory=True)
            run('fit_paper1_solar_lunar_angle.py','--reference-lock',lock,'--north',base/'north',
                '--site',site,'--neutral-separation-deg',70,'--fit-neutral-separation','--output',base/'fit')
        finish(a,run)
        return
    if not all([a.registration_root,a.raw_manifest_root,a.metadata_audit]):p.error('All original reconstruction inputs are required')
    records=json.loads(a.metadata_audit.read_text())
    aliases={'WC_pitch':'WCpitch','ChurcHill':'ChurchHill'}
    selected={}
    handedness={}
    # Preflight both sites before generating any output.
    for site,confirmed_date in SITES.items():
        path=a.registration_root/'runs'/site/'locked_site_references.json'
        lock=json.loads(path.read_text())['sites'][site]
        rows=[r for r in records if aliases.get(Path(r['path']).parts[-2],Path(r['path']).parts[-2])==site]
        if not rows:p.error('No NEF metadata audit records for '+site)
        row=min(rows,key=lambda r:r['metadata']['EXIF DateTimeOriginal'])
        meta=row['metadata']
        moment=datetime.strptime(meta['EXIF DateTimeOriginal'],'%Y:%m:%d %H:%M:%S')
        if moment.date().isoformat()!=confirmed_date:p.error('Confirmed capture date differs: '+site)
        offset=meta.get('EXIF OffsetTimeOriginal')
        if not offset:p.error('Missing explicit NEF UTC offset: '+site)
        observed=datetime.fromisoformat(moment.isoformat()+offset)
        anchor=lock.get('astronomical_anchor') or {}
        if not lock.get('moon_pixel_xy') or not anchor.get('utc'):
            p.error('Reviewed observed Moon and astronomical anchor required: '+site)
        anchored=datetime.fromisoformat(anchor['utc'].replace('Z','+00:00'))
        # Legacy anchors explicitly name this field UTC but omit the suffix.
        if anchored.tzinfo is None:anchored=anchored.replace(tzinfo=timezone.utc)
        if abs((anchored-observed).total_seconds())>1:
            p.error(site+': locked anchor does not match first retained NEF time. '
                    'Correct the metadata-based astronomical anchor and re-lock before running; no silent date substitution.')
        if lock.get('radiometry_status')!='black_subtracted_linear_camera_green':
            p.error('Registered linear-green lock required: '+site)
        selected[site]=path.resolve()
        handedness[site]=lock.get('celestial_orientation',{}).get('handedness','reversed')
    a.output.mkdir(parents=True)
    src=Path(__file__).resolve().parent
    def run(script,*args):
        subprocess.run([sys.executable,str(src/script),*map(str,args)],check=True)
    for site,lock in selected.items():
        base=a.output/site
        base.mkdir()
        (base/'source_lock_path.txt').write_text(str(lock)+'\n')
        run('audit_spider_measured_stokes.py','--reference-lock',lock,'--sites',site,
            '--raw-manifest-root',a.raw_manifest_root,'--output',base/'audit')
        run('rotate_spider_measured_north.py','--reference-lock',lock,
            '--audit',base/'audit'/site,'--site',site,'--handedness',handedness[site],
            '--output',base/'north')
        run('fit_paper1_solar_lunar_angle.py','--reference-lock',lock,
            '--north',base/'north','--site',site,'--neutral-separation-deg','70','--fit-neutral-separation',
            '--output',base/'fit')
    finish(a,run)

    
def finish(a,run):
    run('plot_paper1_two_site_solar_lunar.py','--root',a.output,'--output',a.output/'figure')
    results={site:json.loads((a.output/site/'fit/solar_lunar_fit.json').read_text()) for site in SITES}
    (a.output/'figure/joint_fit_results.json').write_text(json.dumps(results,indent=2)+'\n')
    backup=None
    if a.figure_copy.exists():
        backup=a.figure_copy.with_name(a.figure_copy.name+'_previous_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));a.figure_copy.rename(backup)
    try:archive_figure(a.output/'figure',a.figure_copy)
    except Exception:
        if backup is not None and not a.figure_copy.exists():backup.rename(a.figure_copy)
        raise


if __name__=='__main__':main()

