"""Extend frozen HPS/LED scenarios to historical sites with verified time inputs.

Existing site source scales are preserved. New sites use explicitly conditional
upward-radiance ratios relative to Oxford, not a fitted ground-skyglow relation.
"""
import argparse,csv,hashlib,json,subprocess
from pathlib import Path
from datetime import datetime,timedelta
from zoneinfo import ZoneInfo
import numpy as np
from run_spider_heading_convergence import model_command
ROOT=Path(__file__).resolve().parents[1]


def site_timezone(site_id):
    return 'Europe/Stockholm' if site_id in ['s10','s14'] else ('Europe/Paris' if site_id=='s12' else 'Europe/London')


def extend_scales(config,refs):
    result=json.loads(json.dumps(config));scales=result['alan_scale_by_site']
    oxford=next(r for r in refs if r['site_id']=='s02');oxrad=float(oxford['original_multiyear_monthly_radiance_mean'])
    if oxrad<=0 or 'oxfordbussinesspark' not in scales:raise ValueError('Positive Oxford reference radiance and saved scale required')
    inventory=[]
    for ref in refs:
        site=Path(ref['source_filename']).stem;rad=float(ref['original_multiyear_monthly_radiance_mean'])
        if not np.isfinite(rad) or rad<0:raise ValueError('Invalid site radiance')
        basis='preserved_existing_source_scale' if site in scales else 'conditional_multiyear_upward_radiance_ratio_to_oxford'
        if site not in scales:scales[site]=scales['oxfordbussinesspark']*rad/oxrad
        inventory.append(dict(site=site,site_id=ref['site_id'],radiance_mean=rad,radiance_sd=float(ref['original_multiyear_monthly_radiance_sd']),source_scale=scales[site],scale_basis=basis))
    result['all_site_extension_note']='New site scales use actual reference upward radiance ratios to Oxford, with unchanged HPS/LED DNB factor; conditional scenario, not measured ground skyglow, no DoLP fit.'
    return result,inventory


def prepare_anchors(batch,refs,existing):
    from anchor_moon_from_local_time import capture_time,local_to_utc
    from wang2016_sky import lunar_track
    anchors={};inventory=[]
    for ref in refs:
        site=Path(ref['source_filename']).stem;item=dict(site=site,site_id=ref['site_id'],status='pending_capture_metadata',reason='')
        inventory.append(item)
        if site in existing:
            anchor=dict(existing[site])
            if anchor['site_id']!=ref['site_id']:raise ValueError('Site-ID mismatch in existing anchor')
            # Known inherited anchors are preserved; no new source selection.
            anchors[site]=anchor;item.update(status='ready_existing_ephemeris',timestamp_source=anchor.get('timestamp_source','saved anchor'),timezone=anchor['timezone']);continue
        rp=batch/site/'registration/rotation_analysis.json'
        try:
            report=json.loads(rp.read_text());source=next(r for r in report['inputs'] if Path(r['path']).name=='1.tif')
            stamp,origin=capture_time(Path(source['path']))
            if 'fallback' in origin:raise ValueError('Only TIFF export/fallback time available; independently confirm capture time before use')
            date=datetime.strptime(stamp,'%Y:%m:%d %H:%M:%S').date();reference=datetime.fromisoformat(ref['field_date']).date()
            if date not in [reference,reference+timedelta(days=1)]:raise ValueError('EXIF date does not match field date or following-midnight date')
            zone=site_timezone(ref['site_id']);utc,offset=local_to_utc(stamp,zone)
            el,az=lunar_track([utc],float(ref['latitude']),float(ref['longitude']))
            if not .1<float(el[0])<89.9:raise ValueError('Capture Moon outside supported elevation range; review date/time/location')
            anchors[site]=dict(site_id=ref['site_id'],local_time=stamp,timezone=zone,utc_offset_hours=offset,utc=utc,timestamp_source=origin,
                reference_field_date=ref['field_date'],moon_elevation_deg=float(el[0]),moon_compass_azimuth_deg=float(az[0]),candidate_pixel_xy=None,
                source_report_sha256=hashlib.sha256(rp.read_bytes()).hexdigest(),flags=['Time zone inferred from site country; original camera local-time convention assumed as in the existing campaign','Ephemeris-only forward anchor: no camera bearing selected; not usable for measured north-up overlays'])
            item.update(status='ready_forward_ephemeris_timezone_assumption',timestamp_source=origin,timezone=zone)
        except (FileNotFoundError,KeyError,StopIteration,ValueError) as exc:item['reason']=str(exc)
    return anchors,inventory


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--from-run',type=Path,required=True);p.add_argument('--batch',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--reviewed-anchors',type=Path);p.add_argument('--heading-step',type=int,default=1)
    p.add_argument('--references',type=Path,default=ROOT/'data/reference/prc_spider_18_sites_verified.csv')
    p.add_argument('--prepare-only',action='store_true',help='Write inventory/anchors/scenarios without nightly calculations')
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    old={name:json.loads((a.from_run/name/'provenance.json').read_text()) for name in ['hps','led']}
    existing={r['site']:r['anchor'] for r in old['hps']['sites']}
    if a.reviewed_anchors:
        extra=json.loads(a.reviewed_anchors.read_text())['sites']
        for site,anchor in extra.items():
            if site in existing:
                for key in ['site_id','utc','timezone','moon_elevation_deg','moon_compass_azimuth_deg']:
                    if anchor.get(key)!=existing[site].get(key):p.error('Conflicting ephemeris in supplied anchors: '+site)
            existing[site]=anchor
    refs=list(csv.DictReader(a.references.open()));anchors,status=prepare_anchors(a.batch,refs,existing)
    for name,provenance in old.items():
        for record in provenance['sites']:
            for key in ['utc','timezone','moon_elevation_deg','moon_compass_azimuth_deg']:
                if record['anchor'].get(key)!=anchors[record['site']].get(key):p.error('HPS/LED input ephemerides disagree')
        if not provenance.get('independent_alan_config'):p.error('Saved site ALAN configuration required')
    a.output.mkdir(parents=True)
    anchorpath=a.output/'all_site_ephemeris_anchors.json';anchorpath.write_text(json.dumps(dict(status='FORWARD_EPHEMERIDES_NOT_ALL_CAMERA_ALIGNED',sites=anchors),indent=2)+'\n')
    scales={};inventories={}
    for name,provenance in old.items():
        scales[name],inventories[name]=extend_scales(provenance['independent_alan_config'],refs)
        (a.output/f'{name}_site_alan.json').write_text(json.dumps(scales[name],indent=2)+'\n')
    for row in status:
        for name in old:
            r=next(r for r in inventories[name] if r['site']==row['site'])
            row.update({name+'_source_scale':r['source_scale'],name+'_scale_basis':r['scale_basis']})
        ref=next(r for r in refs if r['site_id']==row['site_id']);row['radiance_mean']=float(ref['original_multiyear_monthly_radiance_mean']);row['radiance_sd']=float(ref['original_multiyear_monthly_radiance_sd']);row['radiance_units']='nW cm^-2 sr^-1'
    with (a.output/'site_input_status.csv').open('w',newline='') as f:
        keys=list(dict.fromkeys(k for r in status for k in r));writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(status)
    (a.output/'extension_notes.json').write_text(json.dumps(dict(ready_sites=list(anchors),pending_sites=[r for r in status if r['site'] not in anchors],
        camera_overlay_status='Additional sites have forward ephemerides; no unverified bright object chosen as Moon',
        scale_note='Actual site radiance used only as explicitly conditional proxy for new source scales; existing Cottage adjustment and Oxford scale preserved',heading_step=a.heading_step,
        reference_sha256=hashlib.sha256(a.references.read_bytes()).hexdigest(),source_provenance_sha256={name:hashlib.sha256((a.from_run/name/'provenance.json').read_bytes()).hexdigest() for name in old}),indent=2)+'\n')
    print(f'{len(anchors)}/{len(refs)} sites have usable forward inputs; inspect {a.output}/site_input_status.csv',flush=True)
    if a.prepare_only:return
    if not anchors:raise ValueError('No usable ephemerides; inventory written')
    for name,provenance in old.items():
        command=model_command(provenance['parameters'],a.output/name,a.heading_step,dict(moon_anchors=str(anchorpath),site_alan_config=str(a.output/f'{name}_site_alan.json'),references=str(a.references),sites=list(anchors)))
        subprocess.run(command,check=True)
    prediction_status=[]
    for name in old:
        provenance=json.loads((a.output/name/'provenance.json').read_text())
        for record in provenance['sites']:
            count=record['heading_records'];prediction_status.append(dict(source=name,site=record['site'],heading_records=count,status='computed' if count else 'no_valid_lunar_field_samples_not_zero_PRC'))
    with (a.output/'site_prediction_status.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(prediction_status[0]));writer.writeheader();writer.writerows(prediction_status)
    levels=old['hps']['parameters']['alan_ratios'];level=5. if 5. in levels else levels[-1]
    fig,ax=plt.subplots(figsize=(max(12,len(anchors)*.8),5),layout='constrained');sites=list(anchors)
    for j,name in enumerate(['hps','led']):
        rows=list(csv.DictReader((a.output/name/'prc_ranges.csv').open()));values=[]
        for site in sites:
            selected=[float(r['median_abs_prc']) for r in rows if r['site']==site and float(r['alan_band_ratio'])==level and float(r['body_tilt_deg'])==0 and float(r['atmospheric_factor'])==old[name]['parameters']['atmospheric_factors'][0]]
            values.append(float(np.median(selected)) if selected else np.nan)
        ax.scatter(np.arange(len(sites))+(j-.5)*.2,values,label=name.upper(),marker=['o','D'][j])
    ax.set_xticks(range(len(sites)),sites,rotation=60,ha='right');ax.set(ylabel='Median over time of heading-median |PRC|',title=f'All ready sites; fixed source scenario level {level:g}; {a.heading_step}deg heading reference\nDifferent capture nights; conditional radiance-based additions at new sites')
    ax.legend();fig.savefig(a.output/'all_sites_prc_overview.png',dpi=170);plt.close(fig)
    print('All ready-site forward scenarios complete; pending metadata remain listed explicitly.',flush=True)

if __name__=='__main__':main()
