"""Summarize visible spectral scenarios and the separate unvalidated UV proxy."""
from pathlib import Path
import argparse,csv,json
from datetime import datetime
import numpy as np
from scipy.integrate import trapezoid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from spider_spectral_lunar import pigment_proxy,ks91_brightness,lunar_phase_distance
from wang2016_sky import berry_normalized_pattern,horizon_factor
from matplotlib.patches import Circle
from scipy.ndimage import map_coordinates
from reconstruct_spider_eye_contrast import digest


def sky_dop_snapshot(provenance,spectra_path,site,atmosphere,alan,grid=241):
    """Capture-time photon-sensitivity-weighted DoLP, not average per-ray DoLP.

    Spectral Q/U integrated before normalization. Uses same relative intensity
    reference and fixed absolute ALAN as nightly model; no fitted AoP offset.
    North up/east right equisolid display. Gray masks KS91-excluded directions.
    """
    record=next(r for r in provenance['sites'] if r['site']==site)
    anchor=record['anchor'];params=provenance['parameters']
    y,x=np.indices((grid,grid),dtype=float);center=(grid-1)/2;radius=center*.96
    dx=(x-center)/radius;dy=-(y-center)/radius;r=np.hypot(dx,dy)
    sky=r<.985;theta=2*np.arcsin(np.clip(r/np.sqrt(2),0,1));phi=np.arctan2(dy,dx)
    el=anchor['moon_elevation_deg'];az=anchor['moon_compass_azimuth_deg'];zs=np.deg2rad(90-el);phim=np.deg2rad(90-az)
    th=theta[sky];ph=phi[sky]
    rho=np.arccos(np.clip(np.cos(th)*np.cos(zs)+np.sin(th)*np.sin(zs)*np.cos(ph-phim),-1,1))
    permitted=rho>=np.deg2rad(10) if params.get('brightness_model','constant')=='ks91' else np.ones(len(th),dtype=bool)
    intensity=np.ones(len(th))
    if params.get('brightness_model','constant')=='ks91':
        phase,distance=lunar_phase_distance([anchor['utc']])
        intensity[:]=np.nan
        intensity[permitted]=ks91_brightness(th[permitted],el,rho[permitted],float(phase[0]),params['extinction_v'],float(distance[0]))/record['capture_zenith_moon_brightness_nL']
    d=np.genfromtxt(spectra_path,delimiter=',',names=True);wave=d['wavelength_nm']
    natural=d['natural']/trapezoid(d['natural'],wave);artificial=d['alan']/trapezoid(d['alan'],wave)
    sensitivity=d['sensitivity'];scale=d['dop_scale'] if 'dop_scale' in d.dtype.names else np.ones(len(wave))
    ni=trapezoid(sensitivity*natural,wave);ai=trapezoid(sensitivity*artificial,wave);nq=trapezoid(sensitivity*natural*scale,wave)
    p=atmosphere*abs(berry_normalized_pattern(th,ph,zs,phim,params['neutral_separation_deg']))*horizon_factor(th)
    noalan=p*nq/ni
    actual_alan=alan*record.get('alan_multiplier_relative_to_capture_zenith',1.)
    mixed=p*intensity*nq/(intensity*ni+actual_alan*ai)
    arrays=[]
    for values in [noalan,mixed]:
        arr=np.full((grid,grid),np.nan);arr[sky]=np.where(permitted,values,np.nan);arrays.append(arr)
    mr=radius*np.sin(zs/2)/np.sin(np.pi/4)
    moon=(center+mr*np.sin(np.deg2rad(az)),center-mr*np.cos(np.deg2rad(az)))
    fov=radius*np.sin(np.deg2rad(15))/np.sin(np.pi/4)
    return arrays,dict(center=center,fov_radius=fov,moon_xy=moon,anchor=anchor,band_nm=[float(wave[0]),float(wave[-1])],actual_alan_capture_zenith_ratio=float(actual_alan))


def measured_north_up(batch,site,geometry,grid=241):
    """Reproject saved experimental Stokes for display using reviewed Moon bearing.

    No registration fit. Camera-channel DoLP remains distinct from proxy-weighted
    model. Bilinear interpolation is for display only; input arrays are untouched.
    """
    folder=batch/site/'registration';rp=folder/'rotation_analysis.json'
    report=json.loads(rp.read_text());anchor=geometry['anchor']
    if digest(rp)!=anchor['source_report_sha256']:raise ValueError('Measured alignment differs from Moon anchor')
    with np.load(folder/'linear_stokes.npz',allow_pickle=False) as z:i,q,u=[z[k] for k in ['I','Q','U']]
    center=np.asarray(report['reference_center_xy_px']);r90=report['radius_60_px']*np.sqrt(2)
    bearing=np.asarray(anchor['candidate_pixel_xy'])-center
    moonphi=np.pi/2-np.deg2rad(anchor['moon_compass_azimuth_deg'])
    rotation=np.arctan2(bearing[1],bearing[0])-moonphi
    y,x=np.indices((grid,grid),dtype=float);c=(grid-1)/2;rad=c*.96
    dx=(x-c)/rad;dy=-(y-c)/rad;r=np.hypot(dx,dy);phi=np.arctan2(dy,dx)+rotation
    coords=[center[1]+r*r90*np.sin(phi),center[0]+r*r90*np.cos(phi)]
    ii,qq,uu=[map_coordinates(a,coords,order=1,mode='constant',cval=np.nan) for a in [i,q,u]]
    with np.errstate(divide='ignore',invalid='ignore'):p=np.hypot(qq,uu)/ii
    physical=(r<.985)&np.isfinite(p)&(ii>0)&(p<=1)
    return np.where(physical,p,np.nan),dict(stokes_sha256=digest(folder/'linear_stokes.npz'),rotation_report_sha256=digest(rp),note='Original measured camera-channel DoLP; spectral response not matched to spider proxy; Moon bearing fixes display orientation')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['spectra','hps','led','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--batch',type=Path,help='Saved experimental all-sites batch; adds measured camera-channel DoLP maps')
    p.add_argument('--atmospheric-factor',type=float,default=.8)
    p.add_argument('--alan-ratio',type=float,default=1.,help='ALAN scenario level; config sets site/spectral scaling')
    a=p.parse_args()
    if not 0<=a.atmospheric_factor<=1 or not np.isfinite(a.alan_ratio) or a.alan_ratio<0:p.error('Invalid atmosphere/ALAN values')
    if a.output.exists():p.error('Use fresh output directory')
    a.output.mkdir(parents=True)
    fig=plt.figure(figsize=(14,16),layout='constrained')
    outer=fig.add_gridspec(2,1,height_ratios=[1.15,1.])
    top=outer[0].subgridspec(2,2)
    axes=np.array([[fig.add_subplot(top[r,c]) for c in range(2)] for r in range(2)])
    ax=axes[0,0]
    for label,name in [('Natural','HPS_visible_input.csv'),('HPS','HPS_visible_input.csv'),('LED 2958–3969K group','LED_mid_visible_input.csv')]:
        d=np.genfromtxt(a.spectra/name,delimiter=',',names=True);wave=d['wavelength_nm']
        values=d['natural'] if label=='Natural' else d['alan']
        ax.plot(wave,values/trapezoid(values,wave),label=label)
    ax.set(title='Source + power-law transfer scenarios; equal band energy',xlabel='Wavelength (nm)',ylabel='Normalized energy spectrum (per nm)');ax.legend(fontsize=8)
    ax=axes[0,1];wave=np.arange(300.,751.)
    for amplitude in [0,.29,1]:ax.plot(wave,pigment_proxy(wave,uv_amplitude=amplitude),label=f'UV amplitude {amplitude:g}')
    ax.axvspan(450,750,color='gray',alpha=.1,label='Current PRC calculation band')
    ax.set(title='510nm pigment proxy; UV strength is a scenario',xlabel='Wavelength (nm)',ylabel='Relative photon sensitivity');ax.legend(fontsize=8)
    for ax,site in zip(axes[1],['Cottage','oxfordbussinesspark']):
        for label,folder in [('HPS',a.hps),('LED 2958–3969K group',a.led)]:
            rows=[r for r in csv.DictReader((folder/'prc_ranges.csv').open()) if r['site']==site and float(r['atmospheric_factor'])==a.atmospheric_factor and float(r['alan_band_ratio'])==a.alan_ratio and float(r['body_tilt_deg'])==0]
            if not rows:raise ValueError('Requested atmosphere/ALAN scenario absent from saved results')
            # Convert UTC to each site's local hour, retaining hours >24.
            from zoneinfo import ZoneInfo
            zone=ZoneInfo('Europe/Stockholm' if site=='Cottage' else 'Europe/London')
            times=[datetime.fromisoformat(r['utc']).astimezone(zone) for r in rows];day=times[0].date()
            hours=[24*(t.date()-day).days+t.hour+t.minute/60 for t in times]
            ax.plot(hours,[float(r['median_abs_prc']) for r in rows],label=label)
        ax.set(title=site,xlabel='Local hour (24 = midnight)',ylabel='Median absolute PRC over body headings',ylim=(0,None));ax.legend(fontsize=8)
    lower=outer[1].subgridspec(2,4 if a.batch else 3)
    skyaxes=[];summary={}
    hp=json.loads((a.hps/'provenance.json').read_text());lp=json.loads((a.led/'provenance.json').read_text())
    for row,site in enumerate(['Cottage','oxfordbussinesspark']):
        hps,geometry=sky_dop_snapshot(hp,a.spectra/'HPS_visible_input.csv',site,a.atmospheric_factor,a.alan_ratio)
        led,other=sky_dop_snapshot(lp,a.spectra/'LED_mid_visible_input.csv',site,a.atmospheric_factor,a.alan_ratio)
        if geometry['anchor']!=other['anchor']:raise ValueError('HPS/LED snapshots use different Moon anchors')
        if hp['parameters'].get('brightness_model')!=lp['parameters'].get('brightness_model') or hp['parameters'].get('extinction_v')!=lp['parameters'].get('extinction_v'):raise ValueError('HPS/LED brightness settings differ')
        summary[site]=dict(geometry,hps_actual_alan_capture_zenith_ratio=geometry['actual_alan_capture_zenith_ratio'],led_actual_alan_capture_zenith_ratio=other['actual_alan_capture_zenith_ratio'])
        panels=[('Natural only',hps[0]),('Natural + HPS',hps[1]),('Natural + LED',led[1])]
        if a.batch:
            observed,obs_meta=measured_north_up(a.batch,site,geometry)
            summary[site]['measured']=obs_meta;panels.insert(0,('Measured camera channel',observed))
        for col,(title,arr) in enumerate(panels):
            ax=fig.add_subplot(lower[row,col]);skyaxes.append(ax);cm=plt.get_cmap('turbo').copy();cm.set_bad('lightgray')
            art=ax.imshow(arr,cmap=cm,vmin=0,vmax=1)
            ax.add_patch(Circle((geometry['center'],geometry['center']),geometry['fov_radius'],fill=False,color='white',lw=1))
            ax.plot(*geometry['moon_xy'],'kx',ms=6)
            ax.set_title(f'{site}: {title}\n{geometry["anchor"]["local_time"]} local',fontsize=10);ax.axis('off')
            if col==0:ax.text(.04,.95,'N ↑',transform=ax.transAxes,ha='left',va='top');ax.text(.96,.5,'E →',transform=ax.transAxes,ha='right',va='center')
    fig.colorbar(art,ax=skyaxes,label='Photon-sensitivity-weighted sky DoLP (450–750nm)',shrink=.75)
    (a.output/'sky_snapshot_settings.json').write_text(json.dumps(dict(atmospheric_factor=a.atmospheric_factor,fixed_alan_capture_zenith_ratio=a.alan_ratio,sites=summary,independent_alan_config=hp.get('independent_alan_config'),note='Natural-only retains atmospheric factor; models integrate proxy-weighted Stokes; measured maps use original camera channel. No DoP fitting'),indent=2)+'\n')
    fig.suptitle('VISIBLE-ONLY pilot: horizontal body, rotating 360deg; KS91 changing moonlight\nAtmospheric factor '+str(a.atmospheric_factor)+'; independent ALAN scenario '+str(a.alan_ratio)+'; sky maps at capture time\nWhite circle = 60deg full field; cross = Moon; measured camera channel differs from model spider weighting; no DoP fit',fontsize=11)
    fig.savefig(a.output/'visible_spectral_nightly_comparison.png',dpi=160);plt.close(fig)

if __name__=='__main__':main()
