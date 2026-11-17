"""Convert audited site RAW sequences in confirmed retained-file acquisition order."""
import argparse,json,re,subprocess,sys
from pathlib import Path
from datetime import datetime
from fractions import Fraction
from prepare_spider_linear_raw import sha256,read_metadata

ALIASES={'ChurcHill':'ChurchHill','WC_pitch':'WCpitch'}
EXCLUSIONS={'museum':{'DSC_5251.NEF':'Recorded 9.4s versus 9.9s; retain prior museum exclusion'},
            'slabhouse':{'DSC_5404.NEF':'Recorded 24.5s versus 24.9s; exclude rather than normalize'}}

def number(value):return float(Fraction(str(value)))
def natural(name):return [int(v) if v.isdigit() else v.casefold() for v in re.split(r'(\d+)',name)]

def make_plans(rows,raw_root):
    grouped={};seen=set()
    for row in rows:
        if row.get('error'):raise ValueError('Audit error: '+row['path'])
        path=Path(row['path']);folder=path.parent.name;filename=path.name
        if path.suffix.lower()!='.nef':raise ValueError('Expected audited NEF input')
        if (folder,filename.casefold()) in seen:raise ValueError('Duplicate audited file')
        seen.add((folder,filename.casefold()));grouped.setdefault(folder,[]).append(row)
    plans=[]
    for folder,rr in sorted(grouped.items()):
        rr=sorted(rr,key=lambda r:natural(Path(r['path']).name));count=len(rr)
        if count not in (4,12):raise ValueError(folder+': expected four or twelve retained files')
        times=[datetime.strptime(r['metadata']['EXIF DateTimeOriginal'],'%Y:%m:%d %H:%M:%S') for r in rr]
        if any(b<=a for a,b in zip(times,times[1:])):raise ValueError(folder+': filename and time order disagree')
        mapping=[];excluded=[];reasons=[];settings=[]
        for i,row in enumerate(rr,1):
            name=Path(row['path']).name;m=row['metadata']
            mapping.append(dict(exposure=i,stem=Path(name).stem,datetime_original=m['EXIF DateTimeOriginal']))
            reason=EXCLUSIONS.get(folder,{}).get(name)
            if reason:excluded.append(i);reasons.append(dict(exposure=i,file=name,reason=reason))
            else:settings.append((m['Image Model'],number(m['EXIF ExposureTime']),number(m.get('EXIF ISOSpeedRatings',m.get('EXIF PhotographicSensitivity'))),number(m['EXIF FNumber'])))
        if len(set(settings))!=1:raise ValueError(folder+': unreviewed exposure/camera settings change')
        if 1 in excluded:raise ValueError('Cannot exclude geometric reference')
        groups=[[1,4,7,10],[2,5,8,11],[3,6,9,12]] if count==12 else [[1,2,3,4]]
        complete=[g for g in groups if not set(g)&set(excluded)]
        if not complete:raise ValueError(folder+': no complete quartet')
        offsets=sorted({r['metadata'].get('EXIF OffsetTimeOriginal','missing') for r in rr})
        if len(offsets)!=1 or offsets[0]=='missing':raise ValueError(folder+': missing or mixed UTC offsets')
        plans.append(dict(site=ALIASES.get(folder,folder),raw_folder=str(raw_root/folder),folder_alias=folder,
            mapping=mapping,expected_exposure_count=count,excluded_exposures=excluded,exclusion_reasons=reasons,
            analyzer_groups=complete,expected_converted_count=count-len(excluded),utc_offset=offsets[0],
            first_local=rr[0]['metadata']['EXIF DateTimeOriginal'],last_local=rr[-1]['metadata']['EXIF DateTimeOriginal'],
            audited_inputs=rr,site_note='Wells Cathedral; separate site, no coordinates inferred' if folder=='wells' else None))
    return plans


def validate_sources(plans):
    # Preflight all sites before decoding any site.
    for plan in plans:
        folder=Path(plan['raw_folder'])
        expected={Path(r['path']).name.casefold() for r in plan['audited_inputs']}
        actual={p.name.casefold() for p in folder.iterdir() if p.is_file() and p.suffix.lower()=='.nef'}
        if actual!=expected:raise ValueError(str(folder)+': inventory changed since metadata audit')
        for row in plan['audited_inputs']:
            path=folder/Path(row['path']).name;m=row['metadata']
            if path.stat().st_size!=row['size_bytes']:raise ValueError('Audited size changed: '+str(path))
            current=read_metadata(path)
            target=dict(datetime_original=m['EXIF DateTimeOriginal'],model=m['Image Model'],
                exposure_seconds=number(m['EXIF ExposureTime']),iso=number(m.get('EXIF ISOSpeedRatings',m.get('EXIF PhotographicSensitivity'))),f_number=number(m['EXIF FNumber']))
            for key,value in target.items():
                if current[key]!=value:raise ValueError('Audited '+key+' changed: '+str(path))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--audit',type=Path,required=True)
    p.add_argument('--raw-root',type=Path,default=Path('data/raw/skyimages'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--channel',choices=['green','blue'],default='green')
    p.add_argument('--sites',nargs='+',help='Optional canonical site names')
    p.add_argument('--plan-only',action='store_true',help='Validate inputs and write mappings/plan without decoding')
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory; existing results are never overwritten')
    try:
        plans=make_plans(json.loads(a.audit.read_text()),a.raw_root)
        if a.sites:
            unknown=set(a.sites)-{r['site'] for r in plans}
            if unknown:raise ValueError('Unknown sites: '+', '.join(sorted(unknown)))
            plans=[r for r in plans if r['site'] in a.sites]
        validate_sources(plans)
    except (ValueError,KeyError,OSError) as e:p.error(str(e))
    a.output.mkdir(parents=True);(a.output/'mappings').mkdir()
    report=dict(status='planned' if a.plan_only else 'running',camera_channel=a.channel,audit_sha256=sha256(a.audit),
        script_sha256=sha256(__file__),sequence_basis='User confirmed retained filename/time order despite deleted filenames',
        sites=[],raw_inputs_modified=False,registration_performed=False)
    for plan in plans:
        site=plan['site'];mapping=a.output/'mappings'/f'{site}.json'
        mapping.write_text(json.dumps(plan['mapping'],indent=2)+'\n')
        result={k:v for k,v in plan.items() if k not in ['audited_inputs','mapping']}
        result['mapping_path']=str(mapping);result['mapping_sha256']=sha256(mapping)
        result['status']='planned'
        if not a.plan_only:
            cmd=[sys.executable,str(Path(__file__).with_name('prepare_spider_linear_raw.py')),
                '--raw-dir',plan['raw_folder'],'--mapping',str(mapping),'--raw-kind','nef',
                '--channel',a.channel,'--output',str(a.output/site)]
            if plan['excluded_exposures']:cmd+=['--exclude-exposures']+list(map(str,plan['excluded_exposures']))
            print('Preparing '+site+'; exclusions '+str(plan['excluded_exposures']),flush=True)
            with (a.output/f'{site}_conversion.log').open('w') as stream:
                process=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT)
            result['returncode']=process.returncode
            result['status']='prepared_not_registered' if process.returncode==0 else 'failed_see_conversion_log'
            if process.returncode==0:
                manifest=a.output/site/'linear_raw_manifest.json';decoded=json.loads(manifest.read_text())
                if len(decoded['frames'])!=plan['expected_converted_count'] or decoded['analyzer_groups']!=plan['analyzer_groups']:
                    result['status']='failed_manifest_disagrees_with_plan'
                result['manifest_path']=str(manifest);result['manifest_sha256']=sha256(manifest)
        report['sites'].append(result)
        (a.output/'batch_summary.json').write_text(json.dumps(report,indent=2)+'\n')
        print(site+': '+result['status'],flush=True)
    failures=[r for r in report['sites'] if r['status'].startswith('failed')]
    report['status']='plan_complete' if a.plan_only else ('partial_failure' if failures else 'conversion_complete_not_registered')
    (a.output/'batch_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Saved '+str(a.output/'batch_summary.json'),flush=True)
    if failures:sys.exit(1)

if __name__=='__main__':main()
