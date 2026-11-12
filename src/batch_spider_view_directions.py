"""Compare zenith, DoLP-band and historical PRC for all 18 spider sites.

Reuses existing Stokes arrays; never performs registration or fits to historic
PRC. Optional per-site run/result overrides support corrected Gavarnie output.
"""
from pathlib import Path
import argparse
import csv
import json
import subprocess
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# Fixed historical plotting order; IDs are checked against the reference CSV at run time.
SITES={
 'ChurchHill':('s07',None),'Cottage':('s14',None),'EastHorrington':('s13',None),
 'Gavarnie':('s12',None),'minifactory':('s05',None),'museum':('s15',None),
 'opensquare':('s16',None),'Ounsdale':('s09',None),'oxfordbussinesspark':('s02',None),
 'PennCC':('s08',None),'Skillinge':('s10',None),'slabhouse':('s17',None),
 'southhorrington':('s11',None),'stannes':('s04',None),'Waitrose':('s03',None),
 'WClawn':('s18',None),'WCpitch':('s01',None),'zoology':('s06',None)}
from reconstruct_spider_eye_contrast import digest
ROOT=Path(__file__).resolve().parents[1]


def overrides(values):
    result={}
    for value in values:
        if '=' not in value: raise ValueError('Use Site=directory for overrides')
        site,directory=value.split('=',1)
        if site not in SITES: raise ValueError(f'Unknown historical site: {site}')
        if site in result: raise ValueError(f'Duplicate override: {site}')
        result[site]=Path(directory).resolve()
    return result


def extract_rows(site,summary,reference,source):
    band=summary['selected_band'];rows=[]
    for r in summary['zenith_and_band_results']:
        row=dict(site=site,site_id=summary['site_id'],method=r['method'],response_mode=r['response_mode'],
            field='zenith' if r['field']=='zenith' else 'band',prc=r['prc'],historical_prc=float(reference),
            difference_from_historical=r['prc']-float(reference),field_elevation_deg=r['elevation_deg'],
            field_azimuth_image_deg=r['azimuth_deg'],mean_dop=r['mean_dop'],physical_fraction=r['physical_fraction'],
            band_selection_elevation_deg=band['elevation_deg'],band_selection_azimuth_image_deg=band['azimuth_deg'],
            full_field_deg=60,angular_radius_deg=30,source_run=str(source),status='REQUIRES_FIELD_AND_ALIGNMENT_REVIEW')
        rows.append(row)
    return rows


def comparison_plot(rows,path):
    order=list(SITES);x=np.arange(len(order))
    fig,axes=plt.subplots(2,1,figsize=(15,9),sharex=True,sharey=True)
    for ax,method in zip(axes,['linear_stokes','legacy_dop_aop']):
        selected={(r['site'],r['field']):r for r in rows if r['method']==method and r['response_mode']=='stokes_intensity'}
        historic={r['site']:r['historical_prc'] for r in rows}
        for offset,field,color,marker,label in [(-.15,'zenith','tab:blue','o','Zenith-centred'),(.15,'band','tab:orange','D','DoLP-band-centred')]:
            values=[selected.get((site,field),{}).get('prc',np.nan) for site in order]
            ax.scatter(x+offset,values,s=42,color=color,marker=marker,label=label,zorder=3)
        ax.scatter(x,[historic.get(site,np.nan) for site in order],s=60,color='black',marker='x',label='Historical',zorder=4)
        ax.set_ylabel('Eye-to-eye log contrast (PRC)');ax.set_title('Direct Stokes' if method=='linear_stokes' else 'Recovered R DoP/AoP averaging')
        ax.grid(axis='y',alpha=.2);ax.legend(ncol=3,fontsize=9)
    axes[-1].set_xticks(x,order,rotation=55,ha='right');axes[-1].set_xlabel('Calibration site')
    fig.suptitle('60deg full field: zenith, DoLP-band centre and historical values\nBand selected by mean DoLP independently of historical PRC; viewing directions require image review')
    fig.tight_layout(rect=[0,0,1,.94])
    for suffix in ['png','svg','pdf']:fig.savefig(path.with_suffix('.'+suffix),dpi=180)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--site-run',action='append',default=[],help='Site=corrected_registration_directory')
    p.add_argument('--reuse-site',action='append',default=[],help='Site=existing_direction_test_directory')
    p.add_argument('--resume',action='store_true')
    a=p.parse_args();batch=a.batch.resolve();out=a.output.resolve()
    if out==batch or batch in out.parents or out in batch.parents:p.error('Output must be separate from batch')
    if out.exists() and not a.resume:p.error('Output exists; use --resume for verified completed site results')
    run_overrides=overrides(a.site_run);reused=overrides(a.reuse_site)
    references={r['site_id']:r for r in csv.DictReader((ROOT/'data/reference/prc_spider_18_sites_verified.csv').open())}
    if set(references)!={v[0] for v in SITES.values()}:p.error('Historical roster and reference IDs disagree')
    out.mkdir(parents=True,exist_ok=True);rows=[];statuses=[]
    for number,(site,(site_id,_)) in enumerate(SITES.items(),1):
        print(f'[{number}/18] {site}: testing field directions',flush=True)
        source=run_overrides.get(site,batch/site/'registration')
        try:
            report=source/'rotation_analysis.json';report_hash=digest(report)
            report_data=json.loads(report.read_text())
            if report_data['site_id']!=site_id:raise ValueError('Source run/site mismatch')
            destination=reused.get(site,out/site)
            summary_path=destination/'view_direction_summary.json'
            if site not in reused and not (a.resume and summary_path.exists()):
                with (out/f'{site}.log').open('w') as log:
                    subprocess.run([sys.executable,str(ROOT/'src/test_spider_view_direction.py'),
                        '--run',str(source),'--output',str(destination)],stdout=log,stderr=subprocess.STDOUT,check=True)
            summary=json.loads(summary_path.read_text())
            if summary['source_report_sha256']!=report_hash or summary['site_id']!=site_id:
                raise ValueError('Existing direction test does not match selected registration run')
            if summary.get('angular_radius_deg')!=30:raise ValueError('Incompatible field radius')
            for method,expected in summary['source_array_sha256'].items():
                if digest(source/f'{method}.npz')!=expected:raise ValueError('Derived Stokes source changed')
            if site not in reused and summary['script_sha256']!=digest(ROOT/'src/test_spider_view_direction.py'):
                raise ValueError('Resume result uses different direction-test code; use fresh output')
            rows.extend(extract_rows(site,summary,references[site_id]['prc_sigg'],source))
            band=summary['selected_band']
            statuses.append(dict(site=site,status='processed_requires_review',result_directory=str(destination),
                band_elevation_deg=band['elevation_deg'],band_azimuth_image_deg=band['azimuth_deg'],
                flagged_exposures=sum(bool(r['flags']) for r in report_data['rotations']),
                review_note='Inspect band overlay for terrain and inherited alignment errors'))
        except (ValueError,OSError,KeyError,subprocess.CalledProcessError) as e:
            statuses.append(dict(site=site,status='failed_or_missing',error=str(e)))
            print(f'  FAILED: {e}',flush=True)
        if rows:
            with (out/'all_sites_view_directions.csv').open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
        (out/'view_direction_batch_status.json').write_text(json.dumps(dict(sites=statuses,
            raw_images_modified=False,historical_values_changed=False,
            band_selection='Maximum mean DoLP; no reference PRC in position selection',
            azimuth='Image-relative; not compass bearing'),indent=2)+'\n')
    if rows:comparison_plot(rows,out/'sites_zenith_band_historical')
    print(f'Finished: comparison plots and CSV in {out}',flush=True)
    if any(s['status']=='failed_or_missing' for s in statuses):sys.exit(1)

if __name__=='__main__':main()
