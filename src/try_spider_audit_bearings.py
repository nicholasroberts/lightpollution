"""Try audit-guided image bearings; preserve all astronomical and PRC inputs."""
import argparse,csv,json,shutil,subprocess,sys
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter
from reconstruct_spider_eye_contrast import digest

ROOT=Path(__file__).resolve().parents[1]
# Approximate directions from the 18-page camera-coordinate audit, 16 Nov 2026.
# Angle from image +x, increasing toward image +y. These are search windows,
# not claimed precise measurements or independent Moon identifications.
SEEDS={'minifactory':-3.,'museum':-174.,'opensquare':-177.,
       'southhorrington':180.,'stannes':-3.,'WCpitch':180.,'zoology':0.}
PRESERVE={'Cottage','oxfordbussinesspark','Skillinge','slabhouse'}


def locate_in_window(intensity,center,r90,elevation,seed,halfwidth=12.,altwidth=6.):
    """Local intensity contrast only; no DoLP/AoP or global source selection."""
    y,x=np.indices(intensity.shape)
    dx=x-center[0];dy=y-center[1];radius=np.hypot(dx,dy)
    angle=np.rad2deg(np.arctan2(dy,dx))
    delta=(angle-seed+180)%360-180
    altitude=90-np.rad2deg(2*np.arcsin(np.clip(radius/(r90*np.sqrt(2)),0,1)))
    valid=np.isfinite(intensity)&(radius<r90*.985)&(intensity>0)
    window=valid&(abs(delta)<=halfwidth)&(abs(altitude-elevation)<=altwidth)
    if not window.any():raise ValueError('No valid pixels in audit-guided search window')
    finite=intensity[valid]
    if not finite.size:raise ValueError('No positive finite sky intensity')
    fill=np.where(valid,intensity,float(np.median(finite)))
    score=gaussian_filter(fill,2)-gaussian_filter(fill,15)
    iy,ix=np.unravel_index(np.argmax(np.where(window,score,-np.inf)),score.shape)
    if not np.isfinite(score[iy,ix]) or score[iy,ix]<=0:
        raise ValueError('No positive local intensity contrast in selected window')
    if abs(delta[iy,ix])>halfwidth-2 or abs(altitude[iy,ix]-elevation)>altwidth-1:
        raise ValueError('Peak touches search-window boundary; leave unaligned for review')
    return [float(ix),float(iy)],dict(
        search_seed_image_angle_deg=seed,search_halfwidth_deg=halfwidth,
        search_elevation_halfwidth_deg=altwidth,selected_image_angle_deg=float(angle[iy,ix]),
        image_implied_elevation_deg=float(altitude[iy,ix]),
        elevation_residual_deg=float(altitude[iy,ix]-elevation),
        method='Audit-guided sector, local intensity contrast; provisional identity, no DoLP/AoP fit')


def window_self_check():
    # A much brighter opposite source must not steal the local selection.
    size=201;center=np.array([100.,100.]);r90=95.;el=30.
    rr=r90*np.sin(np.deg2rad((90-el)/2))/np.sin(np.pi/4)
    y,x=np.indices((size,size));image=np.ones((size,size))
    image+=10*np.exp(-((x-(100+rr))**2+(y-100)**2)/8)
    image+=1000*np.exp(-((x-(100-rr))**2+(y-100)**2)/8)
    xy,meta=locate_in_window(image,center,r90,el,0.)
    if xy[0]<100 or abs(xy[1]-100)>1 or abs(meta['elevation_residual_deg'])>1:
        raise AssertionError('Local source-window check failed')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['run','batch','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use fresh output directory')
    window_self_check()
    old={s:json.loads((a.run/s/'provenance.json').read_text()) for s in ['hps','led']}
    for s,prov in old.items():
        if prov['parameters']['heading_step']!=1:p.error('Requires retained 1 degree results')
    anchors={r['site']:dict(r['anchor']) for r in old['hps']['sites']}
    if set(anchors)!={r['site'] for r in old['led']['sites']}:p.error('Source site lists differ')
    for r in old['led']['sites']:
        if r['anchor']!=anchors[r['site']]:p.error('HPS/LED anchors differ')
    a.output.mkdir(parents=True);display=a.output/'display_run';display.mkdir()
    inventory=[]
    for site,anchor in anchors.items():
        previous=anchor.get('candidate_pixel_xy')
        row=dict(site=site,status='unaligned_pending_review',reason='',previous_pixel_xy=previous)
        # Uncertain sites must not inherit old candidate-based orientation.
        anchor['candidate_pixel_xy']=None
        try:
            if site in PRESERVE or site in SEEDS:
                folder=a.batch/site/'registration';rp=folder/'rotation_analysis.json'
                report=json.loads(rp.read_text())
                if digest(rp)!=anchor.get('source_report_sha256'):
                    raise ValueError('Registration report hash differs from saved anchor')
                center=np.asarray(report['reference_center_xy_px'],dtype=float)
                r90=float(report['radius_60_px'])*np.sqrt(2)
                if site in PRESERVE:
                    if previous is None:raise ValueError('Previously checked bearing missing')
                    xy=np.asarray(previous,dtype=float)
                    if xy.shape!=(2,) or not np.all(np.isfinite(xy)):
                        raise ValueError('Invalid existing bearing')
                    rr=float(np.linalg.norm(xy-center))
                    implied=90-float(np.rad2deg(2*np.arcsin(np.clip(rr/(r90*np.sqrt(2)),0,1))))
                    residual=implied-anchor['moon_elevation_deg']
                    if rr<=0 or rr>=r90 or abs(residual)>6:
                        raise ValueError('Existing candidate radius inconsistent with saved ephemeris')
                    anchor['candidate_pixel_xy']=xy.tolist()
                    row.update(status='preserved_audit_consistent_bearing_provisional',
                               selected_pixel_xy=xy.tolist(),elevation_residual_deg=residual)
                else:
                    with np.load(folder/'linear_stokes.npz',allow_pickle=False) as z:intensity=z['I']
                    xy,meta=locate_in_window(intensity,center,r90,anchor['moon_elevation_deg'],SEEDS[site])
                    anchor['candidate_pixel_xy']=xy
                    row.update(status='audit_guided_local_candidate_provisional',selected_pixel_xy=xy,**meta)
                anchor['display_bearing_review']=dict(status=row['status'],
                    note='Visual audit plus restricted intensity check; requires resulting overlay inspection')
            else:row['reason']='Ambiguous Moon identity; intentionally left in camera coordinates'
        except (OSError,KeyError,ValueError) as exc:
            anchor['candidate_pixel_xy']=None;row['reason']=str(exc)
        inventory.append(row)
        print(site+': '+row['status']+(' — '+row['reason'] if row['reason'] else ''),flush=True)
    for s,prov in old.items():
        derived=json.loads(json.dumps(prov))
        for record in derived['sites']:record['anchor']=anchors[record['site']]
        derived['display_only_alignment']=dict(original_run=str(a.run),
            source_provenance_sha256=digest(a.run/s/'provenance.json'),
            note='Only image-bearing metadata changed. All astronomical values, spectra, ALAN and saved PRCs preserved.')
        target=display/s;target.mkdir()
        (target/'provenance.json').write_text(json.dumps(derived,indent=2)+'\n')
        shutil.copyfile(a.run/s/'prc_ranges.csv',target/'prc_ranges.csv')
        if digest(target/'prc_ranges.csv')!=digest(a.run/s/'prc_ranges.csv'):
            raise AssertionError('PRC results changed during display preparation')
    for name in ['site_input_status.csv','site_prediction_status.csv']:
        if (a.run/name).exists():shutil.copyfile(a.run/name,display/name)
    (a.output/'bearing_review.json').write_text(json.dumps(dict(
        status='PROVISIONAL_DISPLAY_BEARINGS_NOT_ASTRONOMICAL_CORRECTIONS',
        sites=inventory,unaligned_sites=[r['site'] for r in inventory if anchors[r['site']]['candidate_pixel_xy'] is None]),indent=2)+'\n')
    subprocess.run([sys.executable,str(ROOT/'src/audit_spider_moon_positions.py'),
        '--run',str(display),'--batch',str(a.batch),'--output',str(a.output/'bearing_audit')],check=True)
    subprocess.run([sys.executable,str(ROOT/'src/plot_spider_site_details.py'),
        '--run',str(display),'--batch',str(a.batch),'--output',str(a.output/'reports'),
        '--alan-level','5','--atmospheric-factor','0.8','--use-saved-bearings'],check=True)
    print('Inspect bearing_audit/moon_positions_overview.png and reports/all_site_details.pdf',flush=True)


if __name__=='__main__':main()
