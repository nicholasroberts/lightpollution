"""Inspect Moon identity and timing in original camera coordinates; no source fitting."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.backends.backend_pdf import PdfPages
from anchor_moon_from_local_time import capture_time
from reconstruct_spider_eye_contrast import digest


def elevation_radius(r90,elevation):
    return r90*np.sin(np.deg2rad((90-elevation)/2))/np.sin(np.pi/4)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['run','batch','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    prov=json.loads((a.run/'hps/provenance.json').read_text())
    anchors={r['site']:r['anchor'] for r in prov['sites']}
    inventory=a.run/'site_input_status.csv'
    sites=[r['site'] for r in csv.DictReader(inventory.open())] if inventory.exists() else list(anchors)
    a.output.mkdir(parents=True);summary=[];manual={}
    overview,axes=plt.subplots(int(np.ceil(len(sites)/3)),3,figsize=(15,5*int(np.ceil(len(sites)/3))),squeeze=False,layout='constrained')
    for ax in axes.flat:ax.axis('off')
    with PdfPages(a.output/'moon_position_audit.pdf') as pdf:
        for index,site in enumerate(sites):
            anchor=anchors.get(site);item=dict(site=site,saved_anchor=anchor,status='needs_visual_review')
            summary.append(item);manual[site]=dict(confirmed_moon_pixel_xy=None,confirmed_local_time=None,confirmed_timezone=None,notes='')
            fig,axs=plt.subplots(1,3,figsize=(17,7),layout='constrained')
            try:
                folder=a.batch/site/'registration';rp=folder/'rotation_analysis.json'
                report=json.loads(rp.read_text());center=np.array(report['reference_center_xy_px']);r90=report['radius_60_px']*np.sqrt(2)
                with np.load(folder/'linear_stokes.npz',allow_pickle=False) as z:i,q,u=[z[k] for k in ['I','Q','U']]
                yy,xx=np.indices(i.shape);sky=np.hypot(xx-center[0],yy-center[1])<r90
                valid=sky&np.isfinite(i)&(i>0)
                lo,hi=np.percentile(i[valid],[5,99.8]);scaled=np.clip((i-lo)/max(hi-lo,1e-12),0,1)
                with np.errstate(divide='ignore',invalid='ignore'):dop=np.hypot(q,u)/i
                dop=np.where(valid&(dop<=1),dop,np.nan)
                images=[np.where(sky,scaled,np.nan),np.where(sky,np.log1p(100*scaled),np.nan),dop]
                titles=['Intensity: linear stretch','Intensity: logarithmic stretch','Experimental DoLP']
                views=list(zip(axs,images,titles))+[(axes.flat[index],images[1],site)]
                for ax,img,title in views:
                    ax.imshow(img,cmap='turbo' if title=='Experimental DoLP' else 'gray',vmin=0,vmax=1 if title=='Experimental DoLP' else None)
                    if anchor:
                        el=anchor['moon_elevation_deg']
                        ax.add_patch(Circle(center,elevation_radius(r90,el),fill=False,color='cyan',lw=1))
                        candidate=anchor.get('candidate_pixel_xy')
                        if candidate is not None:ax.plot(*candidate,'rx',ms=10,label='Saved candidate; recheck identity')
                    ax.set_title(title);ax.set_xlabel('Original pixel x');ax.set_ylabel('Original pixel y');ax.axis('on')
                raw=next((r for r in report.get('inputs',[]) if Path(r['path']).name=='1.tif'),None)
                try:item['raw_timestamp'],item['raw_timestamp_source']=capture_time(Path(raw['path'])) if raw else (None,'No 1.tif entry')
                except (OSError,ValueError) as exc:item['raw_timestamp_error']=str(exc)
                item['report_sha256']=digest(rp)
                if anchor:
                    item['report_matches_anchor']=anchor.get('source_report_sha256')==item['report_sha256']
                    candidate=anchor.get('candidate_pixel_xy')
                    if candidate is not None:
                        rr=np.linalg.norm(np.asarray(candidate)-center)
                        item['candidate_implied_elevation_deg']=float(90-np.rad2deg(2*np.arcsin(np.clip(rr/(r90*np.sqrt(2)),0,1))))
                    text=(f"Saved local time: {anchor['local_time']} | zone: {anchor['timezone']} | UTC: {anchor['utc']}\n"
                          f"Calculated Moon elevation {anchor['moon_elevation_deg']:.2f}°; azimuth {anchor['moon_compass_azimuth_deg']:.2f}°\n"
                          f"Original metadata: {item.get('raw_timestamp','unavailable')} ({item.get('raw_timestamp_source',item.get('raw_timestamp_error',''))})")
                else:text='No usable saved ephemeris: inspect capture date/time and Moon identity before modelling.'
                fig.suptitle(site+' — original camera coordinates; camera north not assumed\n'+text+'\nCyan: predicted elevation ring. Red: saved image candidate, NOT verified by this audit. No source selected or fitted.',fontsize=11)
            except (OSError,KeyError,ValueError) as exc:
                item.update(status='missing_or_invalid_inputs',reason=str(exc));fig.suptitle(site+': '+str(exc))
                for ax in axs:ax.axis('off')
                axes.flat[index].text(.5,.5,site+'\nInputs unavailable',ha='center',va='center')
            fig.savefig(a.output/f'{site}_moon_audit.png',dpi=140);pdf.savefig(fig);plt.close(fig)
    overview.suptitle('Moon-position audit in camera coordinates\nCyan: predicted elevation ring; red: saved candidate requiring visual confirmation',fontsize=15)
    overview.savefig(a.output/'moon_positions_overview.png',dpi=140);plt.close(overview)
    (a.output/'moon_audit_details.json').write_text(json.dumps(summary,indent=2)+'\n')
    (a.output/'moon_review_template.json').write_text(json.dumps(manual,indent=2)+'\n')
    print('Saved camera-coordinate audit; no Moon positions or model parameters changed.',flush=True)


if __name__=='__main__':main()
