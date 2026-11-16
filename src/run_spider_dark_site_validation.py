"""Apply recorded visible Moon bearings; rerun all four sites with frozen AoP offsets.

Raw Stokes, ephemeris elevations and animal-model polarization basis are unchanged.
"""
import argparse,csv,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def corrected_anchors(sources,overrides,centers=None):
    merged={}
    for source in sources:
        for site,anchor in source['sites'].items():
            if site in merged and merged[site]!=anchor:raise ValueError(f'Conflicting anchors: {site}')
            merged[site]=dict(anchor)
    for site,override in overrides['sites'].items():
        if site not in merged:raise ValueError(f'Missing ephemeris for {site}')
        point=np.asarray(override['candidate_pixel_xy'],dtype=float)
        if point.shape!=(2,) or not np.all(np.isfinite(point)):raise ValueError('Finite x,y coordinates required')
        anchor=merged[site]
        if override.get('bearing_rotation_deg',0):
            if centers is None or site not in centers:raise ValueError('Saved camera centre required for bearing rotation')
            point=rotate_point(point,centers[site],override['bearing_rotation_deg'])
        anchor['previous_candidate_pixel_xy']=anchor['candidate_pixel_xy']
        anchor['candidate_pixel_xy']=point.tolist()
        anchor['bearing_override']=dict(override,source=overrides['source'],status=overrides['status'])
        anchor['flags']=[f for f in anchor.get('flags',[]) if not f.startswith('Ring-constrained brightness candidate:')]
        anchor['flags'].append('Approximate visible-source bearing from uploaded intensity overlay; centroid and camera tilt remain uncertain')
    return dict(sites=merged,status='EPHEMERIS_ELEVATION_WITH_EXPLICIT_VISUAL_BEARING_OVERRIDES')


def rotate_point(point,center,degrees):
    angle=np.deg2rad(float(degrees));c=np.cos(angle);s=np.sin(angle)
    return np.asarray(center)+np.array([[c,-s],[s,c]])@(np.asarray(point)-np.asarray(center))


def profile_fwhm(x,y):
    """Zero-baseline half maximum around main peak; no extrapolation/smoothing."""
    x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float)
    if len(x)<3 or x.shape!=y.shape or np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)) or np.any(np.diff(x)<=0) or np.any(y<0) or np.max(y)<=0:
        return dict(status='invalid_or_insufficient_profile',fwhm_deg=None)
    peak=int(np.argmax(y));half=y[peak]/2;left=right=None
    for j in range(peak-1,-1,-1):
        if y[j]<=half:
            if x[j+1]-x[j]<=15:left=float(x[j]+(half-y[j])*(x[j+1]-x[j])/(y[j+1]-y[j]))
            break
    for j in range(peak+1,len(x)):
        if y[j]<=half:
            if x[j]-x[j-1]<=15:right=float(x[j-1]+(half-y[j-1])*(x[j]-x[j-1])/(y[j]-y[j-1]))
            break
    complete=left is not None and right is not None
    return dict(status='complete' if complete else 'crossing_not_resolved_in_sampled_support',
        peak_deg=float(x[peak]),peak_value=float(y[peak]),half_maximum=float(half),
        left_crossing_deg=left,right_crossing_deg=right,
        left_hwhm_deg=float(x[peak]-left) if left is not None else None,
        right_hwhm_deg=float(right-x[peak]) if right is not None else None,
        fwhm_deg=float(right-left) if complete else None,
        sampled_support_deg=[float(x[0]),float(x[-1])],
        note='Binned means; linear crossing interpolation, no baseline subtraction, no confidence interval; absent crossings are not extrapolated')


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from compare_wang2016_measured_skies import ephemeris_source_pixel
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',type=Path,required=True)
    p.add_argument('--moon-anchors',type=Path,nargs='+',required=True)
    p.add_argument('--fixed-offset-audit',type=Path,required=True)
    p.add_argument('--bearing-overrides',type=Path,default=ROOT/'config/spider_dark_site_moon_bearings_20261113.json')
    p.add_argument('--neutral-separation-deg',type=float,default=10.)
    p.add_argument('--sky-masks',type=Path)
    p.add_argument('--capture-alan-ratios',type=Path)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    frozen=json.loads(a.fixed_offset_audit.read_text())
    if 'neutral_separation_deg' in frozen:
        if frozen['neutral_separation_deg']!=a.neutral_separation_deg:p.error('Neutral separation differs from frozen audit')
    sources=[json.loads(f.read_text()) for f in a.moon_anchors]
    centers={site:json.loads((a.batch/site/'registration/rotation_analysis.json').read_text())['reference_center_xy_px'] for source in sources for site in source['sites']}
    payload=corrected_anchors(sources,json.loads(a.bearing_overrides.read_text()),centers)
    required=['Cottage','oxfordbussinesspark','Skillinge','slabhouse']
    if any(site not in payload['sites'] for site in required):p.error('All four sites require saved ephemeris anchors')
    # Verify provenance and pixel bounds BEFORE creating output or running audit.
    reports={}
    for site in required:
        path=a.batch/site/'registration/rotation_analysis.json';raw=path.read_bytes();report=json.loads(raw)
        if hashlib.sha256(raw).hexdigest()!=payload['sites'][site]['source_report_sha256']:p.error(f'Anchor alignment hash mismatch: {site}')
        with np.load(a.batch/site/'registration/linear_stokes.npz',allow_pickle=False) as z:shape=z['I'].shape
        x,y=payload['sites'][site]['candidate_pixel_xy']
        if not 0<=x<shape[1] or not 0<=y<shape[0]:p.error(f'Candidate outside saved image: {site}')
        reports[site]=report
    a.output.mkdir(parents=True)
    anchorpath=a.output/'corrected_moon_anchors.json';anchorpath.write_text(json.dumps(payload,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(13,6),layout='constrained')
    for ax,site in zip(axes,['Skillinge','slabhouse']):
        anchor=payload['sites'][site];report=reports[site];center=report['reference_center_xy_px'];r90=report['radius_60_px']*np.sqrt(2)
        with np.load(a.batch/site/'registration/linear_stokes.npz',allow_pickle=False) as z:i=z['I']
        y,x=np.indices(i.shape);sky=np.hypot(x-center[0],y-center[1])<r90
        ax.imshow(np.where(sky,i,np.nan),cmap='gray',vmin=np.nanpercentile(i[sky],5),vmax=np.nanpercentile(i[sky],99.9))
        ax.plot(*anchor['candidate_pixel_xy'],'r+',ms=12,label='Visible-source bearing estimate')
        projected=ephemeris_source_pixel(anchor,center,r90)
        ax.plot(*projected,'cx',ms=10,label='Model source at ephemeris elevation')
        radius=np.linalg.norm(np.array(projected)-center)
        ax.add_patch(plt.Circle(center,radius,fill=False,color='cyan'))
        ax.plot(*anchor['previous_candidate_pixel_xy'],'yo',mfc='none',label='Rejected automatic edge')
        ax.set_title(f'{site}: Moon elevation {anchor["moon_elevation_deg"]:.2f} deg');ax.legend(fontsize=8)
    fig.suptitle('Corrected visible-source bearings; radial position fixed by ephemeris\nApproximate overlay estimates, not fitted to DoLP or AoP')
    fig.savefig(a.output/'corrected_moon_review.png',dpi=160);plt.close(fig)
    cmd=[sys.executable,str(ROOT/'src/audit_spider_aop_conventions.py'),'--batch',str(a.batch),'--moon-anchors',str(anchorpath),'--reviewed-sites',*required,'--fixed-offset-audit',str(a.fixed_offset_audit),'--neutral-separation-deg',str(a.neutral_separation_deg),'--output',str(a.output/'audit')]
    for flag,value in [('--sky-masks',a.sky_masks),('--capture-alan-ratios',a.capture_alan_ratios)]:
        if value:cmd.extend([flag,str(value)])
    subprocess.run(cmd,check=True)
    inventory={r['site']:r for r in csv.DictReader((a.output/'audit/site_radiance_inventory.csv').open())}
    missing=[site for site in required if inventory[site]['status']!='diagnostic_ready']
    if missing:raise RuntimeError(f'Audit incomplete for {missing}; inspect inventory; no complete summary claimed')
    profile=list(csv.DictReader((a.output/'audit/dark_site_dop_profiles.csv').open()))
    fig,axes=plt.subplots(1,3,figsize=(15,4),layout='constrained')
    ratios=json.loads(a.capture_alan_ratios.read_text()) if a.capture_alan_ratios else {}
    for ax,site in zip(axes,['Cottage','Skillinge','slabhouse']):
        rows=[r for r in profile if r['site']==site];angle=np.array([float(r['moon_separation_deg']) for r in rows]);natural=np.array([float(r['natural_eta1_dop']) for r in rows])
        ax.plot(angle,[float(r['measured_dop']) for r in rows],color='black',label='Measured camera')
        for eta in [1.,.8,.6,.4]:ax.plot(angle,natural*eta/(1+ratios.get(site,0)),label=f'eta={eta}')
        ax.set(xlabel='Moon separation (deg)',ylabel='Mean DoLP',ylim=(0,1),title=f'{site}; radiance {float(inventory[site]["radiance_mean"]):.3f}')
        ax.legend(fontsize=8)
    fig.suptitle('Three comparison sites: prescribed atmospheric factors; original measured radiances retained\n'+('Conditional supplied camera-band ALAN ratios' if ratios else 'Natural-only curves; ALAN unspecified, not assumed absent from measurements'))
    fig.savefig(a.output/'three_site_dop_comparison.png',dpi=170);plt.close(fig)
    widths={};fig,axes=plt.subplots(1,3,figsize=(15,4),layout='constrained')
    for ax,site in zip(axes,['Cottage','Skillinge','slabhouse']):
        rows=[r for r in profile if r['site']==site]
        angle=np.array([float(r['moon_separation_deg']) for r in rows]);measured=np.array([float(r['measured_dop']) for r in rows]);natural=np.array([float(r['natural_eta1_dop']) for r in rows])
        widths[site]={'measured':profile_fwhm(angle,measured),'natural':profile_fwhm(angle,natural)}
        for values,label,color in [(measured,'Measured','black'),(natural,'Natural model','tab:blue')]:
            if len(values) and max(values)>0:ax.plot(angle,values/max(values),label=label,color=color)
        ax.axhline(.5,color='gray',ls='--');ax.set(xlabel='Moon separation (deg)',ylabel='Peak-normalized DoLP',ylim=(0,1.05),title=site);ax.legend()
    fig.suptitle('Band width independent of uniform amplitude scaling; dashed = half maximum')
    fig.savefig(a.output/'three_site_normalized_profiles.png',dpi=170);plt.close(fig)
    (a.output/'profile_widths.json').write_text(json.dumps(widths,indent=2,allow_nan=False)+'\n')
    audit=json.loads((a.output/'audit/convention_and_atmosphere_audit.json').read_text())
    summary=dict(status='FOUR_SITE_CAMERA_CONVENTION_AND_DOP_AUDIT',selected_convention='conjugate',fixed_offset_deg=audit['common_offsets_deg']['conjugate'],
        aop={s:audit['aop_metrics'][s]['conjugate'] for s in required},radiances={s:float(inventory[s]['radiance_mean']) for s in required},
        source_bearing_uncertainty='Skillinge/slabhouse approximate +/-25 raw pixels from visual overlay; review corrected_moon_review.png',
        profile_widths=widths,physical_model_basis_changed=False,atmospheric_factors_fitted=False,alan_inferred_from_radiance=False,
        input_sha256={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in [*a.moon_anchors,a.fixed_offset_audit,a.bearing_overrides]},command=cmd)
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print('Complete:',a.output/'summary.json')

if __name__=='__main__':main()
