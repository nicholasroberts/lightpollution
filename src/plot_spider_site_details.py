"""Per-site spectral, nightly PRC and capture-time DoLP reports; no refitting."""
import argparse,csv,json
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np
from scipy.integrate import trapezoid
from scipy.ndimage import map_coordinates
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from plot_spider_spectral_pilot import sky_dop_snapshot,measured_north_up,fov_statistics
from reconstruct_spider_eye_contrast import digest
from spider_plot_time import date_coordinates,format_time_axis


def measured_camera_frame(batch,site,grid=241):
    """Equal-solid-angle camera-frame display when north has not been reviewed."""
    folder=batch/site/'registration'
    report=json.loads((folder/'rotation_analysis.json').read_text())
    with np.load(folder/'linear_stokes.npz',allow_pickle=False) as z:
        arrays=[z[k] for k in ['I','Q','U']]
    center=report['reference_center_xy_px'];r90=report['radius_60_px']*np.sqrt(2)
    y,x=np.indices((grid,grid),dtype=float);c=(grid-1)/2;rad=c*.96
    dx=(x-c)/rad;dy=(y-c)/rad
    coords=[center[1]+dy*r90,center[0]+dx*r90]
    i,q,u=[map_coordinates(v,coords,order=1,mode='constant',cval=np.nan) for v in arrays]
    with np.errstate(divide='ignore',invalid='ignore'):p=np.hypot(q,u)/i
    return np.where((np.hypot(dx,dy)<.985)&(i>0)&(p<=1),p,np.nan),dict(
        orientation='native camera frame; north unknown; not spatially aligned to model',
        stokes_sha256=digest(folder/'linear_stokes.npz'))


def selected_rows(rows,site,atmosphere,level):
    return sorted([r for r in rows if r['site']==site and
        float(r['atmospheric_factor'])==atmosphere and
        float(r['alan_band_ratio'])==level and float(r['body_tilt_deg'])==0],key=lambda r:r['utc'])


def locked_measured(lockpath,site,anchor,grid=241):
    """Read frozen products; never rerun registration or infer north here."""
    root=lockpath.parent;manifest=json.loads(lockpath.read_text())
    checks=json.loads((root/'reference_sha256.json').read_text())
    if digest(lockpath)!=checks.get(lockpath.name):raise ValueError('Reference manifest hash mismatch')
    record=manifest['sites'][site]
    if record['status']!='numerically_locked_awaiting_visual_review':
        raise ValueError('Registration unresolved: '+record['status'])
    for key in ['utc','timezone','moon_elevation_deg','moon_compass_azimuth_deg']:
        if (record.get('astronomical_anchor') or {}).get(key)!=anchor.get(key):
            raise ValueError('Locked astronomical anchor differs: '+key)
    center=np.asarray(record['reference_center_xy_px']);radius=record['rim_radius_px']
    north=record.get('north_up_dop_path')
    if north:
        path=root/north
        if digest(path)!=checks.get(str(path.relative_to(root))):raise ValueError('North display hash mismatch')
        with np.load(path,allow_pickle=False) as z:values=z['DoLP']
        orientation='North up; Moon-based pure rotation; provisional identity'
    else:
        path=root/site/'registration/linear_stokes.npz'
        if digest(path)!=checks.get(str(path.relative_to(root))):raise ValueError('Stokes product hash mismatch')
        with np.load(path,allow_pickle=False) as z:
            with np.errstate(divide='ignore',invalid='ignore'):values=np.hypot(z['Q'],z['U'])/z['I']
        orientation='Registered camera frame; north unresolved'
    y,x=np.indices((grid,grid),dtype=float);c=(grid-1)/2;r=c*.96
    dx=(x-c)/r;dy=(y-c)/r
    arr=map_coordinates(values,[center[1]+dy*radius,center[0]+dx*radius],order=1,mode='constant',cval=np.nan)
    arr=np.where((np.hypot(dx,dy)<.985)&(arr>=0)&(arr<=1),arr,np.nan)
    return arr,dict(reference_lock_sha256=digest(lockpath),orientation=orientation,
        registration_status=record['status'],moon_localization=record.get('moon_localization'),
        note='Frozen products; numerical landmark checks passed; visual review still required'),bool(north)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True,help='all_sites directory containing hps/ and led/')
    p.add_argument('--batch',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--alan-level',type=float,default=5.)
    p.add_argument('--atmospheric-factor',type=float,default=.8)
    p.add_argument('--use-saved-bearings',action='store_true',help='Opt in to saved candidate image bearings; these are not independently verified by this script')
    p.add_argument('--reference-lock',type=Path,help='Frozen locked_site_references.json; use these registered experimental products and fixed north rotations')
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    if a.reference_lock:
        if not a.reference_lock.is_file():
            p.error('Reference lock not found: '+str(a.reference_lock)+'; run lock_spider_site_registration.py first or supply its actual output path. No figures generated.')
        if not (a.reference_lock.parent/'reference_sha256.json').is_file():
            p.error('Reference checksum file missing beside '+str(a.reference_lock)+'; registration output is incomplete. No figures generated.')
    sources=['hps','led'];provs={s:json.loads((a.run/s/'provenance.json').read_text()) for s in sources}
    spectra={};rows={}
    for s,prov in provs.items():
        par=prov['parameters']
        if par['heading_step']!=1:p.error('Reports require the retained 1 degree heading grid')
        if a.alan_level not in par['alan_ratios'] or a.atmospheric_factor not in par['atmospheric_factors']:p.error('Requested scenario absent')
        path=Path(par['spectra']) if par.get('spectra') else None
        if path is None:p.error('Supplied spectral inputs required')
        if digest(path)!=prov['spectra']['sha256']:p.error('Spectral input changed since calculation')
        spectra[s]=(path,np.genfromtxt(path,delimiter=',',names=True))
        rows[s]=list(csv.DictReader((a.run/s/'prc_ranges.csv').open()))
    for key in ['natural_dop_law','neutral_separation_deg','brightness_model','extinction_v']:
        if provs['hps']['parameters'].get(key)!=provs['led']['parameters'].get(key):p.error('Source physical settings differ: '+key)
    sites=[r['site'] for r in provs['hps']['sites']]
    if set(sites)!={r['site'] for r in provs['led']['sites']}:p.error('Source site lists differ')
    a.output.mkdir(parents=True);summaries=[]
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(a.output/'all_site_details.pdf') as pdf:
        for site in sites:
            records={s:next(r for r in provs[s]['sites'] if r['site']==site) for s in sources}
            if records['hps']['anchor']!=records['led']['anchor']:raise ValueError('Different source capture anchors: '+site)
            anchor=records['hps']['anchor'];zone=ZoneInfo(anchor['timezone'])
            capture=datetime.fromisoformat(anchor['utc'])
            if capture.tzinfo is None:capture=capture.replace(tzinfo=ZoneInfo('UTC'))
            capture=capture.astimezone(zone)
            fig=plt.figure(figsize=(14,12),layout='constrained')
            gs=fig.add_gridspec(4,4,height_ratios=[1,1.05,1.35,.5])
            ax=fig.add_subplot(gs[0,:2])
            for label,s,col in [('Moonlight','hps','natural'),('HPS','hps','alan'),('LED','led','alan')]:
                d=spectra[s][1];wave=d['wavelength_nm'];v=d[col]
                ax.plot(wave,v/trapezoid(v,wave),label=label)
            ax.set(title='Spectral shapes after transfer; equal band energy',xlabel='Wavelength (nm)',ylabel='Normalized energy per nm');ax.legend()
            ax=fig.add_subplot(gs[0,2:])
            for s in sources:
                d=spectra[s][1];ax.plot(d['wavelength_nm'],d['sensitivity'],label=s.upper()+' input sensitivity')
            ax.set(title='510 nm green-sensitive pigment; no added UV peak',xlabel='Wavelength (nm)',ylabel='Relative photon sensitivity');ax.legend()
            ax=fig.add_subplot(gs[1,:])
            for s,color in [('hps','tab:blue'),('led','tab:orange')]:
                selected=selected_rows(rows[s],site,a.atmospheric_factor,a.alan_level)
                if selected:
                    times=[]
                    for r in selected:
                        t=datetime.fromisoformat(r['utc'])
                        if t.tzinfo is None:t=t.replace(tzinfo=ZoneInfo('UTC'))
                        times.append(t.astimezone(zone))
                    ax.plot(date_coordinates(times),[float(r['median_abs_prc']) for r in selected],color=color,label=s.upper())
            ax.axvline(date_coordinates([capture])[0],color='gray',ls='--',label='Image capture')
            ax.set(title=f'ALAN level {a.alan_level:g}; 1° headings; horizontal body; 60° full field',ylabel='Heading-median |PRC|',ylim=(0,None))
            format_time_axis(ax,zone);ax.legend(ncol=3)
            snapshots={};geometry=None
            for s in sources:
                snapshots[s],g=sky_dop_snapshot(provs[s],spectra[s][0],site,a.atmospheric_factor,a.alan_level)
                if geometry is None:geometry=g
            aligned=a.use_saved_bearings and anchor.get('candidate_pixel_xy') is not None
            observed=None;obsmeta={};obsstatus='unavailable'
            try:
                if a.reference_lock:
                    observed,obsmeta,aligned=locked_measured(a.reference_lock,site,anchor)
                    obsstatus=obsmeta['orientation']
                else:
                    if aligned:observed,obsmeta=measured_north_up(a.batch,site,geometry)
                    else:observed,obsmeta=measured_camera_frame(a.batch,site)
                    obsstatus='saved bearing; orientation provisional' if aligned else 'camera frame; north unknown'
            except (FileNotFoundError,ValueError,KeyError) as exc:
                obsmeta={'reason':str(exc)};obsstatus='unavailable: '+str(exc)
            panels=[('Measured camera DoLP',observed),('Natural only',snapshots['hps'][0]),('Natural + HPS',snapshots['hps'][1]),('Natural + LED',snapshots['led'][1])]
            stats={};mapaxes=[];cm=plt.get_cmap('turbo').copy();cm.set_bad('lightgray')
            for j,(title,arr) in enumerate(panels):
                ax=fig.add_subplot(gs[2,j]);mapaxes.append(ax)
                if arr is None:
                    ax.text(.5,.5,'Measured Stokes unavailable\nSee report status',ha='center',va='center',transform=ax.transAxes);ax.axis('off');continue
                art=ax.imshow(arr,cmap=cm,vmin=0,vmax=1)
                ax.add_patch(Circle((geometry['center'],geometry['center']),geometry['fov_radius'],fill=False,color='white'))
                if j>0 or aligned:
                    ax.plot(*geometry['moon_xy'],'kx');ax.text(.03,.97,'N ↑',transform=ax.transAxes,va='top')
                stats[title]=fov_statistics(arr,geometry)
                st=stats[title];label='No valid FoV pixels' if st['mean'] is None else f"FoV mean {st['mean']:.3f}; median {st['median']:.3f}"
                ax.text(.5,-.04,label,transform=ax.transAxes,ha='center',va='top',fontsize=9)
                ax.set_title(title+('\n'+obsstatus if j==0 else '\nNorth up'),fontsize=10);ax.axis('off')
            fig.colorbar(art,ax=mapaxes,shrink=.8,label='DoLP (0–1); camera channel and model weighting differ')
            ref=records['hps']['reference'];rad=float(ref['original_multiyear_monthly_radiance_mean'])
            ratios={s:a.alan_level*records[s].get('alan_multiplier_relative_to_capture_zenith',1.) for s in sources}
            scale_basis='Conditional radiance-ratio source scaling; existing Cottage/Oxford settings retained'
            info=(f"Capture: {capture.strftime('%d %b %Y %H:%M:%S %Z')} | Site radiance: {rad:.4g} nW cm⁻² sr⁻¹ (historical upward measurement)\n"
                  f"Moon elevation {anchor['moon_elevation_deg']:.2f}°; azimuth {anchor['moon_compass_azimuth_deg']:.2f}° | Atmospheric factor {a.atmospheric_factor:g}\n"
                  f"ALAN / natural band energy at capture zenith: HPS {ratios['hps']:.3g}; LED {ratios['led']:.3g} | {scale_basis}\n"
                  'White circle: 60° full FoV. Maps at image time. Measured camera spectral response differs from model pigment; no DoLP fit.')
            ax=fig.add_subplot(gs[3,:]);ax.axis('off');ax.text(0,1,info,va='top',fontsize=9)
            fig.suptitle(site+' — spectra, nightly PRC and capture-time DoLP\n'+provs['hps']['parameters']['natural_dop_law']+f'; ALAN scenario {a.alan_level:g}',fontsize=14)
            fig.savefig(a.output/f'{site}_spectral_prc_dop.png',dpi=160);pdf.savefig(fig);plt.close(fig)
            summaries.append(dict(site=site,capture_local=capture.isoformat(),radiance=rad,measured_status=obsstatus,measured_metadata=obsmeta,fov_statistics=stats,alan_capture_ratios=ratios))
            print('Saved site report: '+site,flush=True)
    (a.output/'site_report_details.json').write_text(json.dumps(dict(heading_step=1,alan_level=a.alan_level,atmospheric_factor=a.atmospheric_factor,sites=summaries),indent=2)+'\n')


if __name__=='__main__':main()
