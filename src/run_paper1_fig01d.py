"""Run existing audited Stokes/north/angle-fit stages for Figure 1D; no registration refit."""
import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

SITES={'PennCC':'2019-04-20','WCpitch':'2019-05-14'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--registration-root',type=Path,required=True)
    p.add_argument('--raw-manifest-root',type=Path,required=True)
    p.add_argument('--metadata-audit',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    records=json.loads(a.metadata_audit.read_text())
    aliases={'WC_pitch':'WCpitch','ChurcHill':'ChurchHill'}
    selected={}
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
            '--audit',base/'audit'/site,'--site',site,'--handedness','reversed',
            '--output',base/'north')
        run('fit_paper1_solar_lunar_angle.py','--reference-lock',lock,
            '--north',base/'north','--site',site,'--neutral-separation-deg','24',
            '--output',base/'fit')
    run('plot_paper1_two_site_solar_lunar.py','--root',a.output,'--output',a.output/'figure')


if __name__=='__main__':main()
