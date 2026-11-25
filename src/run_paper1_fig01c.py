"""Joint shared-separation/solar-fraction fits for Cottage and museum, initial82deg."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from datetime import datetime,timezone
from paper1_figure_archive import archive_figure


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cottage-root',type=Path,default=Path('data/processed/paper1_cottage_north_01'))
    p.add_argument('--museum-north',type=Path,default=Path('data/processed/paper1_museum_north_01'))
    p.add_argument('--museum-lock',type=Path,default=Path('data/processed/spider_registration_green_all_sites_01/runs/museum/locked_site_references.json'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--figure-copy',type=Path,required=True)
    p.add_argument('--replace-figure-copy',action='store_true')
    a=p.parse_args()
    if a.output.exists() or (a.figure_copy.exists() and not a.replace_figure_copy):p.error('Use fresh output; replacing figure-copy requires --replace-figure-copy')
    src=Path(__file__).resolve().parent
    def run(script,*args):subprocess.run([sys.executable,str(src/script),*map(str,args)],check=True)
    sites=[('Cottage',a.cottage_root/'registration/locked_site_references.json',a.cottage_root/'north'),('museum',a.museum_lock,a.museum_north)]
    a.output.mkdir(parents=True);panels=a.output/'panels';panels.mkdir();summaries={}
    for site,lock,north in sites:
        fit=a.output/site/'fit'
        run('fit_paper1_solar_lunar_angle.py','--reference-lock',lock,'--north',north,'--site',site,
            '--neutral-separation-deg',82,'--fit-neutral-separation','--output',fit)
        run('plot_paper1_solar_lunar_panel.py','--fit',fit,'--reference-lock',lock,'--north',north,'--output',panels/site)
        summaries[site]=json.loads((fit/'solar_lunar_fit.json').read_text())
    (panels/'joint_fit_results.json').write_text(json.dumps(summaries,indent=2)+'\n')
    backup=None
    if a.figure_copy.exists():
        backup=a.figure_copy.with_name(a.figure_copy.name+'_previous_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));a.figure_copy.rename(backup)
    try:archive_figure(panels,a.figure_copy)
    except Exception:
        if backup is not None and not a.figure_copy.exists():backup.rename(a.figure_copy)
        raise
    print('Saved both site comparisons:',a.figure_copy)


if __name__=='__main__':main()
