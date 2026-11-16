"""Heading-grid check from frozen previous-run settings, then optional all-site run.

Evaluate 1deg once and exactly subsample the nested 5/10deg grids. At a given
heading, field quadrature is independent of heading-grid spacing. No smoothing.
"""
import argparse,csv,hashlib,json,subprocess,sys
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def model_command(parameters,output,heading_step=1,overrides=None):
    values=dict(parameters);values.update(overrides or {});values.update(output=str(output),heading_step=heading_step)
    cmd=[sys.executable,str(ROOT/'src/model_spider_nightly_prc.py')]
    names=['brightness_model','extinction_v','moon_anchors','output','prc_colour_limit','site_alan_config','spectra','references','atmospheric_factors','alan_ratios','tilts','step_minutes','heading_step','samples','neutral_separation_deg','natural_dop_law','sites']
    for key in names:
        v=values.get(key)
        if v is None:continue
        cmd.append('--'+key.replace('_','-'));cmd.extend(map(str,v if isinstance(v,list) else [v]))
    return cmd


def grid_statistics(headings,values,steps=(10,5,1)):
    headings=np.asarray(headings,dtype=int);values=np.asarray(values,dtype=float)
    if len(headings)!=360 or set(headings)!=set(range(360)) or np.any(~np.isfinite(values)):raise ValueError('Complete finite 0..359deg heading grid required')
    result={}
    for step in steps:
        selected=abs(values[headings%step==0])
        result[step]=dict(min_abs_prc=float(min(selected)),median_abs_prc=float(np.median(selected)),max_abs_prc=float(max(selected)))
    return result


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from spider_plot_time import date_coordinates,format_time_axis
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--from-run',type=Path,required=True,help='Existing normalized run with hps/led provenance')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--add-sites',action='store_true')
    p.add_argument('--batch',type=Path)
    p.add_argument('--reviewed-anchors',type=Path,help='Corrected four-site anchors, optional for all-site metadata reuse')
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    if a.add_sites and not a.batch:p.error('--add-sites requires --batch')
    if a.reviewed_anchors and not a.reviewed_anchors.is_file():p.error('Reviewed anchor file missing')
    if a.batch and not a.batch.is_dir():p.error('Measured batch directory missing')
    saved={name:json.loads((a.from_run/name/'provenance.json').read_text()) for name in ['hps','led']}
    for name,old in saved.items():
        if old['parameters'].get('natural_dop_law')!='rayleigh_normalized':p.error('Use a normalized source run')
        params=old['parameters']
        if params.get('site_alan_config'):
            if json.loads(Path(params['site_alan_config']).read_text())!=old.get('independent_alan_config'):p.error('Saved ALAN input has changed: '+name)
        if params.get('spectra') and old.get('spectra',{}).get('sha256'):
            if hashlib.sha256(Path(params['spectra']).read_bytes()).hexdigest()!=old['spectra']['sha256']:p.error('Saved spectrum has changed: '+name)
    a.output.mkdir(parents=True)
    all_stats=[];commands=[];summary={};groups_by_source={}
    for source,old in saved.items():
        folder=a.output/source
        cmd=model_command(old['parameters'],folder,1);commands.append(cmd)
        print(f'{source}: evaluating 1deg headings; time/field/ALAN settings unchanged',flush=True)
        subprocess.run(cmd,check=True)
        groups=defaultdict(list)
        for row in csv.DictReader((folder/'heading_responses.csv').open()):
            key=(row['site'],row['utc'],row['body_tilt_deg'],row['atmospheric_factor'],row['alan_band_ratio'])
            groups[key].append((int(row['body_heading_deg']),float(row['signed_prc'])))
        if not groups:raise ValueError('No valid lunar samples to compare')
        groups_by_source[source]={}
        for key,values in groups.items():
            statistics=grid_statistics(*zip(*values));groups_by_source[source][key]=statistics
            for step,stat in statistics.items():
                row=dict(source=source,site=key[0],utc=key[1],body_tilt_deg=key[2],atmospheric_factor=key[3],alan_band_ratio=key[4],heading_step_deg=step,**stat)
                for measure in stat:row['difference_vs_1deg_'+measure]=stat[measure]-statistics[1][measure]
                all_stats.append(row)
        errors={}
        for step in [10,5]:
            errors[str(step)]={}
            for measure in ['min_abs_prc','median_abs_prc','max_abs_prc']:
                values=[abs(r['difference_vs_1deg_'+measure]) for r in all_stats if r['source']==source and r['heading_step_deg']==step]
                errors[str(step)][measure]=dict(max_abs_difference=float(max(values)),p95_abs_difference=float(np.quantile(values,.95)))
        summary[source]=errors
        # Same scenarios/times, only heading grid changes. One panel per site.
        sites=list(dict.fromkeys(k[0] for k in groups));fig,axes=plt.subplots(1,len(sites),figsize=(7*len(sites),4.5),squeeze=False,layout='constrained')
        params=old['parameters'];level=5. if 5. in params['alan_ratios'] else params['alan_ratios'][-1]
        tilt=str(float(params['tilts'][0]));eta=str(float(params['atmospheric_factors'][0]))
        for ax,site in zip(axes.flat,sites):
            keys=sorted([k for k in groups if k[0]==site and float(k[2])==float(tilt) and float(k[3])==float(eta) and float(k[4])==level],key=lambda k:k[1])
            stamps=[datetime.fromisoformat(k[1]) for k in keys];x=date_coordinates(stamps)
            for step in [10,5,1]:ax.plot(x,[groups_by_source[source][k][step]['median_abs_prc'] for k in keys],label=f'{step}deg headings')
            zone=ZoneInfo(next(r['anchor']['timezone'] for r in old['sites'] if r['site']==site));format_time_axis(ax,zone)
            ax.set(title=f'{site}; {source.upper()}; ALAN level {level:g}',ylabel='Median absolute PRC');ax.legend()
        fig.suptitle('Heading-grid check: identical time/field/spectrum/ALAN; no smoothing')
        fig.savefig(a.output/f'{source}_heading_convergence.png',dpi=170);plt.close(fig)
    with (a.output/'heading_grid_statistics.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(all_stats[0]));w.writeheader();w.writerows(all_stats)
    (a.output/'convergence_summary.json').write_text(json.dumps(dict(status='GRID_SENSITIVITY_VERSUS_1DEG_NOT_PROOF_OF_CONTINUOUS_CONVERGENCE',errors=summary,commands=commands,
        reference='1deg grid; 5/10deg are exact nested subsamples of the same per-heading calculations',
        statement='Keep 1deg for the extension; no physical parameter changes or data smoothing. Field/time convergence is a separate question.',
        input_sha256={name:hashlib.sha256((a.from_run/name/'provenance.json').read_bytes()).hexdigest() for name in saved}),indent=2)+'\n')
    print('Heading-grid check complete:',a.output/'convergence_summary.json',flush=True)
    if a.add_sites:
        command=[sys.executable,str(ROOT/'src/extend_spider_nightly_sites.py'),'--from-run',str(a.from_run),'--batch',str(a.batch),'--heading-step','1','--output',str(a.output/'all_sites')]
        if a.reviewed_anchors:command+=['--reviewed-anchors',str(a.reviewed_anchors)]
        subprocess.run(command,check=True)

if __name__=='__main__':main()
