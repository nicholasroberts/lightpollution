"""Create consistently scaled per-site PRC figures and an all-site status table."""
import argparse
import csv
import json
from pathlib import Path
from types import SimpleNamespace
from plot_spider_prc_equivalence_figure import run


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True)
    p.add_argument('--reference-lock',type=Path,required=True)
    p.add_argument('--inventory-lock',type=Path,help='Optional full-site registration lock for pending-site accounting')
    p.add_argument('--fit-root',type=Path,required=True)
    p.add_argument('--prc-root',type=Path,required=True)
    p.add_argument('--camera-response',type=Path,default=Path(__file__).resolve().parents[1]/'data/reference/d850_elf_response.csv')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use fresh output directory')
    available={f.parent.name for f in a.results.glob('*/prc_lamp_equivalence.json')}
    locked=json.loads(a.reference_lock.read_text())['sites']
    inventory=json.loads(a.inventory_lock.read_text())['sites'] if a.inventory_lock else locked
    sites=sorted(set(inventory)|set(locked)|available)
    reports=[json.loads((a.results/site/'prc_lamp_equivalence.json').read_text()) for site in available]
    positive=[r['results']['RH1_510']['fits'][k]['energy'] for r in reports for k in ['HPS','LED'] if r['results']['RH1_510']['fits'][k]['energy'] is not None]
    energy_limit=max(positive)*1.2 if positive else 1.
    prc_limit=max([.7]+[r['results']['RH1_510']['fits']['HPS']['natural_prc']*1.05 for r in reports])
    a.output.mkdir(parents=True);status=[]
    for site in sites:
        if site not in available:
            status.append(dict(site=site,status='pending_upstream_PRC_equivalence',detail='No completed spectral PRC-equivalence report; registration/Moon/RAW requirements must be resolved upstream.'))
            continue
        try:
            run(SimpleNamespace(results=a.results,output=a.output/site,sites=[site],camera_response=a.camera_response,
                reference_lock=a.reference_lock,fit_root=a.fit_root,prc_root=a.prc_root,energy_limit=energy_limit,prc_limit=prc_limit))
            status.append(dict(site=site,status='figure_generated',detail='Existing fitted results and locked inputs verified; no new fitting or registration approval.'))
        except (OSError,ValueError,KeyError) as exc:
            status.append(dict(site=site,status='failed_validation',detail=str(exc)))
    with (a.output/'site_figure_status.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['site','status','detail']);writer.writeheader();writer.writerows(status)
    (a.output/'site_figure_status.json').write_text(json.dumps(status,indent=2)+'\n')
    print(json.dumps(status,indent=2))
    if any(v['status']=='failed_validation' for v in status):raise SystemExit(2)


if __name__=='__main__':main()
