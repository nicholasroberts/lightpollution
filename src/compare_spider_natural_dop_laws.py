"""Three-site matched-pixel DoLP-law comparison; no atmosphere/ALAN fitting.

Candidate W/(1+cos²gamma) retains Berry phase/zeros but is a hypothesis, not a
re-derivation of Wang eq17 or a calibrated multiple-scattering Stokes solution.
"""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from wang2016_sky import berry_normalized_pattern,horizon_factor
from run_spider_dark_site_validation import profile_fwhm
ROOT=Path(__file__).resolve().parents[1]
SITES=('Cottage','Skillinge','slabhouse')


def rayleigh_normalized_candidate(w,gamma):
    return np.asarray(w)/(1+np.cos(np.asarray(gamma))**2)


def matched_profiles(gamma,observed,current,candidate,valid):
    rows=[]
    for low in range(10,180,10):
        take=valid&(gamma>=max(15,low))&(gamma<low+10)
        if take.sum()<20:continue
        # Keep the previous fixed bin centres and Moon exclusion: 15..135 etc.
        rows.append(dict(moon_separation_deg=low+5,pixels=int(take.sum()),
            measured_dop=float(np.mean(observed[take])),current_dop=float(np.mean(current[take])),candidate_dop=float(np.mean(candidate[take]))))
    return rows


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',type=Path,required=True)
    p.add_argument('--moon-anchors',type=Path,required=True,help='Corrected anchors from spider_dark_site_validation_02')
    p.add_argument('--sky-masks',type=Path,help='Same SITE.npy boolean sky masks as prior audit, if supplied')
    p.add_argument('--neutral-separation-deg',type=float,default=10.)
    p.add_argument('--reference',type=Path,default=ROOT/'data/reference/prc_spider_18_sites_verified.csv')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    if not 0<a.neutral_separation_deg<90:p.error('Separation must be between 0 and 90deg')
    anchors=json.loads(a.moon_anchors.read_text())['sites']
    references={Path(r['source_filename']).stem:r for r in csv.DictReader(a.reference.open())}
    loaded={};summary={};allrows=[]
    # Validate every site's inputs before writing an output directory.
    for site in SITES:
        folder=a.batch/site/'registration';rp=folder/'rotation_analysis.json';sp=folder/'linear_stokes.npz'
        anchor=anchors[site];report=json.loads(rp.read_text())
        if sha(rp)!=anchor['source_report_sha256']:p.error(f'Moon anchor belongs to another registration: {site}')
        if not 0<anchor['moon_elevation_deg']<90:p.error(f'Moon elevation outside supported range: {site}')
        with np.load(sp,allow_pickle=False) as z:i,q,u=[z[k].astype(float) for k in ['I','Q','U']]
        if i.ndim!=2 or q.shape!=i.shape or u.shape!=i.shape:p.error('Stokes array shape mismatch')
        center=np.asarray(report['reference_center_xy_px']);r90=float(report['radius_60_px'])*np.sqrt(2)
        source=np.asarray(anchor['candidate_pixel_xy'])-center
        if not np.all(np.isfinite(source)) or np.linalg.norm(source)==0:p.error('Invalid source bearing')
        az=np.arctan2(source[1],source[0]);zs=np.deg2rad(90-anchor['moon_elevation_deg'])
        y,x=np.indices(i.shape);r=np.hypot(x-center[0],y-center[1])/r90
        theta=2*np.arcsin(np.clip(r/np.sqrt(2),0,1));phi=np.arctan2(y-center[1],x-center[0])
        gamma=np.arccos(np.clip(np.cos(theta)*np.cos(zs)+np.sin(theta)*np.sin(zs)*np.cos(phi-az),-1,1))
        domain=(r<.985)&(theta<=np.deg2rad(60))
        maskpath=a.sky_masks/f'{site}.npy' if a.sky_masks else None
        if a.sky_masks and not maskpath.exists():p.error(f'Requested sky mask missing: {site}')
        if maskpath:
            mask=np.load(maskpath,allow_pickle=False)
            if mask.dtype!=bool or mask.shape!=i.shape:p.error(f'Mask shape/type mismatch: {site}')
            domain &=mask
        w=np.full(i.shape,np.nan+1j*np.nan,dtype=complex)
        w[domain]=berry_normalized_pattern(theta[domain],phi[domain],zs,az,a.neutral_separation_deg)*horizon_factor(theta[domain])
        candidate=abs(rayleigh_normalized_candidate(w,gamma));current=abs(w)
        with np.errstate(divide='ignore',invalid='ignore'):measured=np.hypot(q,u)/i
        physical=np.isfinite(measured)&(i>0)&(measured<=1)&np.isfinite(current)&(current<=1)&np.isfinite(candidate)&(candidate<=1)
        valid=domain&physical&(gamma>=np.deg2rad(15))
        if valid.sum()<20:p.error(f'Too few matched physical pixels: {site}')
        rows=matched_profiles(np.rad2deg(gamma),measured,current,candidate,valid)
        if len(rows)<3:p.error(f'Insufficient angular support: {site}')
        for row in rows:allrows.append(dict(site=site,**row))
        angle=np.array([r['moon_separation_deg'] for r in rows]);profiles={key:np.array([r[key] for r in rows]) for key in ['measured_dop','current_dop','candidate_dop']}
        if any(max(v)<=0 for v in profiles.values()):p.error(f'Zero-peak profile cannot be normalized: {site}')
        normalized={k:v/max(v) for k,v in profiles.items()}
        metrics={key:dict(width=profile_fwhm(angle,values),
            peak_mean_dop=float(max(values)),
            profile_shape_rmse=float(np.sqrt(np.mean((normalized[key]-normalized['measured_dop'])**2)))) for key,values in profiles.items()}
        cap=valid&(theta<=np.deg2rad(30))
        summary[site]=dict(radiance_mean=float(references[site]['original_multiyear_monthly_radiance_mean']),
            radiance_sd=float(references[site]['original_multiyear_monthly_radiance_sd']),radiance_units='nW cm^-2 sr^-1',
            matched_pixels=int(valid.sum()),excluded_nonphysical_pixels=int((domain&~physical).sum()),
            zenith_fov_30deg_radius={key:dict(mean=float(np.mean(v[cap])),median=float(np.median(v[cap]))) for key,v in [('measured',measured),('current',current),('candidate',candidate)]} if cap.any() else None,
            profiles=metrics,mask_status='external_sky_mask' if maskpath else 'NO_TERRAIN_MASK_provisional',
            mask_sha256=sha(maskpath) if maskpath else None,stokes_sha256=sha(sp),registration_sha256=sha(rp),anchor=anchor)
        loaded[site]=dict(rows=rows,profiles=profiles,normalized=normalized,angle=angle,domain=domain,valid=valid,measured=measured,current=current,candidate=candidate,center=center,r90=r90)
    a.output.mkdir(parents=True)
    fig,axes=plt.subplots(2,3,figsize=(15,8),layout='constrained')
    for col,site in enumerate(SITES):
        d=loaded[site]
        for key,label,color in [('measured_dop','Measured camera','black'),('current_dop','Current Berry/horizon envelope','tab:blue'),('candidate_dop','Rayleigh-denominator candidate','tab:orange')]:
            axes[0,col].plot(d['angle'],d['profiles'][key],label=label,color=color)
            axes[1,col].plot(d['angle'],d['normalized'][key],label=label,color=color)
        axes[0,col].set(title=f'{site}; radiance {summary[site]["radiance_mean"]:.3f}',ylabel='Mean DoLP',ylim=(0,1))
        axes[1,col].axhline(.5,color='gray',ls='--');axes[1,col].set(xlabel='Moon separation (deg)',ylabel='Peak-normalized DoLP',ylim=(0,1.05))
        axes[0,col].legend(fontsize=8)
    fig.suptitle('Matched pixels and bins: two natural DoLP laws versus measurements; no atmosphere or ALAN fit\nCandidate preserves Berry phase/zeros; finite-splitting extension is a hypothesis, not full Wang eq17')
    fig.savefig(a.output/'three_site_natural_law_comparison.png',dpi=170);plt.close(fig)
    for site,d in loaded.items():
        fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
        for ax,key,label in zip(axes,['measured','current','candidate'],['Measured camera','Current envelope','Normalized candidate']):
            im=ax.imshow(np.where(d['valid'],d[key],np.nan),vmin=0,vmax=1,cmap='turbo')
            ax.add_patch(plt.Circle(d['center'],d['r90']*np.sin(np.deg2rad(15))/np.sin(np.pi/4),fill=False,color='white'))
            ax.set_title(label);ax.set_axis_off()
        fig.colorbar(im,ax=list(axes),label='DoLP',shrink=.7);fig.suptitle(f'{site}: identical compared pixels; white circle = 60deg full field')
        fig.savefig(a.output/f'{site}_matched_dop_maps.png',dpi=150);plt.close(fig)
    with (a.output/'matched_profiles.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(allrows[0]));writer.writeheader();writer.writerows(allrows)
    result=dict(status='CANDIDATE_COMPARISON_NOT_PRODUCTION_CHANGE',sites=summary,neutral_separation_deg=a.neutral_separation_deg,
        atmospheric_factor=1.,alan_added=0.,production_model_changed=False,fitted_parameters=[],
        anchor_sha256=sha(a.moon_anchors),reference_sha256=sha(a.reference),script_sha256=sha(Path(__file__)),
        assumptions=['Candidate complex fractional polarization = current W H / (1+cos²gamma); phase and zeros unchanged.',
        'This is an explicit finite-neutral-splitting extension motivated by ideal Rayleigh normalization; not a derived multiple-scattering solution.',
        'Same physical pixels, mask, 15deg Moon exclusion, 60deg-radius audit domain and 10deg bins for all laws; previous bin centres retained.',
        'Atmospheric factor is fixed at 1 and no ALAN added to theoretical curves; real measurements are not assumed ALAN-free.',
        'Peak normalization is only a separate shape diagnostic. Absolute predictions and zenith-field statistics are unscaled.',
        'Shape RMSE uses equal angular-bin weights, not independent-pixel uncertainty. FWHM uses linear bin interpolation without extrapolation.',
        'Historical camera spectral channel differs from spider green510 weighting. Site radiance does not directly specify capture-time ground ALAN.',
        'No AoP offset estimated or applied: phase is identical between candidates, amplitude comparison needs no AoP calibration.'])
    (a.output/'comparison.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Comparison saved:',a.output)

if __name__=='__main__':main()
