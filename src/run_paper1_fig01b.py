"""Recalculate the two-site lunar baseline from accepted RAW celestial locks."""
import argparse
from datetime import datetime, timezone
import json
import subprocess
import sys
from pathlib import Path
from paper1_figure_archive import archive_figure, digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cottage-root',type=Path,default=Path('data/processed/paper1_cottage_north_01'))
    p.add_argument('--museum-north',type=Path,default=Path('data/processed/paper1_museum_north_01'))
    p.add_argument('--museum-lock',type=Path,default=Path('data/processed/spider_registration_green_all_sites_01/runs/museum/locked_site_references.json'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--figure-copy',type=Path,required=True)
    p.add_argument('--replace-figure-copy',action='store_true',help='Retain existing figure copy in a timestamped backup, then replace it')
    a=p.parse_args()
    if a.output.exists() or (a.figure_copy.exists() and not a.replace_figure_copy):
        p.error('Output and figure-copy must be fresh versioned directories')
    cottage_lock=a.cottage_root/'registration/locked_site_references.json'
    anchor=json.loads(cottage_lock.read_text())['sites']['Cottage']['astronomical_anchor']
    if anchor['utc']!='2018-08-21T20:23:49Z':
        p.error('Expected verified Cottage RAW reference UTC 2018-08-21T20:23:49Z')
    sites=[('Cottage',cottage_lock,a.cottage_root/'north'),('museum',a.museum_lock,a.museum_north)]
    # Verify both input lineages before starting expensive model calculations.
    for site,lock,north in sites:
        nr=json.loads((north/'north_rotation.json').read_text())
        if nr['site']!=site or nr['source_lock_sha256']!=digest(lock):
            p.error(site+': north/reference-lock mismatch')
    a.output.mkdir(parents=True)
    src=Path(__file__).resolve().parent
    def run(script,*args):
        subprocess.run([sys.executable,str(src/script),*map(str,args)],check=True)
    for site,lock,north in sites:
        run('model_paper1_lunar_baseline.py','--reference-lock',lock,'--north',north,
            '--site',site,'--neutral-separation-deg',82,'--output',a.output/site)
    run('plot_paper1_lunar_two_sites.py','--cottage',a.output/'Cottage',
        '--cottage-lock',cottage_lock,'--museum',a.output/'museum',
        '--museum-lock',a.museum_lock,'--output',a.output/'panel')
    # Include complete model reports alongside the publication artwork/stats.
    reports={site:json.loads((a.output/site/'lunar_only.json').read_text()) for site,_,_ in sites}
    (a.output/'panel/model_provenance.json').write_text(json.dumps(dict(
        sites=reports,runner_sha256=digest(__file__),neutral_separation_deg=82,
        model='Existing Berry/Wang kernel; no ALAN, fitted dilution or DoLP ceiling',
        source_model_arrays={site:dict(path=str(a.output/site/'lunar_only_stokes.npz'),
            sha256=digest(a.output/site/'lunar_only_stokes.npz')) for site,_,_ in sites}
    ),indent=2)+'\n')
    backup=None
    if a.figure_copy.exists():
        stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        backup=a.figure_copy.with_name(a.figure_copy.name+'_before_82deg_'+stamp)
        a.figure_copy.rename(backup)
    try:
        archive_figure(a.output/'panel',a.figure_copy)
    except Exception:
        if backup is not None and not a.figure_copy.exists():backup.rename(a.figure_copy)
        raise
    if backup is not None:print('Previous figure snapshot retained:',backup)
    print('Saved recalculated 1B and verified snapshot',a.figure_copy)


if __name__=='__main__':main()
