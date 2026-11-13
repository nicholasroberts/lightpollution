"""Cross-site AoP convention audit; retain independent site radiances for DoLP.

Cottage alone estimates one common angle offset per declared basis hypothesis.
Other sites test those fixed offsets. This diagnoses conventions; it does not
establish an optical calibration or infer atmosphere from DoLP alone.
"""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DARK_SITES=('Cottage','Skillinge','slabhouse')
HYPOTHESES=('direct','conjugate','direct_radial_plus','direct_radial_minus',
            'conjugate_radial_plus','conjugate_radial_minus')


def transform(z,phi,name):
    """Explicit spin-2 alternatives, without moving any measured image pixels."""
    if name not in HYPOTHESES:raise ValueError('Unknown convention')
    value=np.conj(z) if name.startswith('conjugate') else z
    if name.endswith('_plus'):value=value*np.exp(2j*phi)
    if name.endswith('_minus'):value=value*np.exp(-2j*phi)
    return value


def fit_offset(obs,model):
    return .5*np.angle(np.mean(obs/abs(obs)*np.conj(model/abs(model))))


def residual_stats(obs,model):
    delta=np.rad2deg(.5*np.angle(obs*np.conj(model)))
    unit=obs/abs(obs)*np.conj(model/abs(model))
    return dict(pixels=len(obs),mae_deg=float(np.mean(abs(delta))),
        median_abs_deg=float(np.median(abs(delta))),p90_abs_deg=float(np.quantile(abs(delta),.9)),
        signed_mean_axial_deg=float(np.rad2deg(.5*np.angle(np.mean(unit)))),
        axial_coherence=float(abs(np.mean(unit))))


def apparent_eta(observed,natural,alan_ratio=None):
    """Eta only conditional on independently supplied capture ALAN/natural ratio.

    Missing ALAN leaves an apparent attenuation, NOT an atmospheric estimate.
    Uniform ratio is a deliberately simple diagnostic, not a radiance transfer.
    """
    if alan_ratio is not None and (not np.isfinite(alan_ratio) or alan_ratio<0):raise ValueError('ALAN ratio must be finite and nonnegative')
    ratio=observed/natural
    result=dict(observed_over_natural_quantiles=np.quantile(ratio,[.1,.5,.9]).tolist())
    if alan_ratio is not None:
        eta=ratio*(1+alan_ratio)
        result.update(conditional_eta_quantiles=np.quantile(eta,[.1,.5,.9]).tolist(),
            capture_alan_over_natural=alan_ratio,fraction_conditional_eta_above_one=float(np.mean(eta>1)))
    else:result['interpretation']='Atmosphere and ALAN remain confounded; not an atmospheric factor estimate'
    return result


def main():
    from compare_wang2016_measured_skies import compare_site
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',type=Path,required=True)
    p.add_argument('--moon-anchors',type=Path,nargs='+',required=True)
    p.add_argument('--reviewed-sites',nargs='+',default=['Cottage','oxfordbussinesspark'],help='Only independently reviewed Moon candidates are admitted')
    p.add_argument('--calibration-site',default='Cottage')
    p.add_argument('--neutral-separation-deg',type=float,required=True)
    p.add_argument('--sky-masks',type=Path)
    p.add_argument('--capture-alan-ratios',type=Path,help='JSON {SITE: capture ALAN/natural ratio}, independently justified, same camera band; never inferred from VIIRS directly')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--reference',type=Path,default=ROOT/'data/reference/prc_spider_18_sites_verified.csv')
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    anchors={}
    for path in a.moon_anchors:
        for site,anchor in json.loads(path.read_text())['sites'].items():
            if site in anchors and anchors[site]!=anchor:p.error(f'Conflicting Moon anchors for {site}')
            anchors[site]=anchor
    ratios=json.loads(a.capture_alan_ratios.read_text()) if a.capture_alan_ratios else {}
    for value in ratios.values():
        if not isinstance(value,(int,float)) or not np.isfinite(value) or value<0:p.error('Invalid ALAN ratio')
    refs={Path(r['source_filename']).stem:r for r in csv.DictReader(a.reference.open())}
    loaded={};inventory=[]
    for site in dict.fromkeys([*DARK_SITES,*a.reviewed_sites]):
        ref=refs[site]
        item=dict(site=site,dark_comparison_site=site in DARK_SITES,
            radiance_mean=float(ref['original_multiyear_monthly_radiance_mean']),
            radiance_sd=float(ref['original_multiyear_monthly_radiance_sd']),
            radiance_ratio_to_cottage=float(ref['original_multiyear_monthly_radiance_mean'])/float(refs['Cottage']['original_multiyear_monthly_radiance_mean']),
            radiance_units='nW cm^-2 sr^-1',status='pending_independent_moon_review')
        inventory.append(item)
        if site not in anchors or site not in a.reviewed_sites:continue
        try:
            baseline,arrays=compare_site(a.batch/site/'registration',None,a.neutral_separation_deg,anchor=anchors[site],alan_relative=0.)
            intensity,obs,model,sky,field,center,r90=arrays
            model=model*np.exp(-2j*np.deg2rad(baseline['global_aop_offset_deg']))
            y,x=np.indices(obs.shape);phi=np.arctan2(y-center[1],x-center[0])
            theta=2*np.arcsin(np.clip(np.hypot(x-center[0],y-center[1])/(r90*np.sqrt(2)),0,1))
            domain=sky&(theta<=np.deg2rad(60))
            maskpath=a.sky_masks/f'{site}.npy' if a.sky_masks else None
            if maskpath and maskpath.exists():
                m=np.load(maskpath,allow_pickle=False)
                if m.dtype!=bool or m.shape!=obs.shape:raise ValueError('Sky mask shape/type mismatch')
                domain &= m;item['mask_status']='external_mask'
            else:item['mask_status']='no_terrain_mask_provisional'
            physical=np.isfinite(obs)&np.isfinite(model)&(intensity>0)&(abs(obs)<=1)
            moonphi=np.arctan2(*(np.asarray(anchors[site]['candidate_pixel_xy'])-center)[::-1])
            zs=np.deg2rad(90-anchors[site]['moon_elevation_deg'])
            gamma=np.rad2deg(np.arccos(np.clip(np.cos(theta)*np.cos(zs)+np.sin(theta)*np.sin(zs)*np.cos(phi-moonphi),-1,1)))
            # Fixed 15 degree exclusion avoids the direct lunar/glare source.
            valid=domain&physical&(gamma>=15)
            ao=valid&(abs(obs)>=.1)&(abs(model)>=.1)
            if ao.sum()<20:raise ValueError('Too few polarized pixels')
            loaded[site]=dict(obs=obs,model=model,phi=phi,ao=ao,valid=valid,field=field,center=center,r90=r90,gamma=gamma)
            item['status']='diagnostic_ready';item['stokes_sha256']=baseline['stokes_sha256']
            item['registration_flagged_exposures']=baseline['flagged_exposures']
            item['sky_mask_sha256']=hashlib.sha256(maskpath.read_bytes()).hexdigest() if maskpath and maskpath.exists() else ''
        except (FileNotFoundError,ValueError,KeyError) as exc:item.update(status='blocked',reason=str(exc))
    if a.calibration_site not in loaded:p.error('Calibration site missing or not independently reviewed')
    a.output.mkdir(parents=True)
    cal=loaded[a.calibration_site];sel=cal['ao']
    offsets={name:float(fit_offset(cal['obs'][sel],transform(cal['model'][sel],cal['phi'][sel],name))) for name in HYPOTHESES}
    metrics={};dop={};profile_rows=[]
    for site,d in loaded.items():
        metrics[site]={};fig,axes=plt.subplots(2,3,figsize=(12,8),layout='constrained')
        for ax,name in zip(axes.flat,HYPOTHESES):
            prediction=transform(d['model'],d['phi'],name)*np.exp(2j*offsets[name])
            metrics[site][name]=residual_stats(d['obs'][d['ao']],prediction[d['ao']])
            residual=np.rad2deg(.5*np.angle(d['obs']*np.conj(prediction)))
            im=ax.imshow(np.where(d['ao'],residual,np.nan),cmap='twilight',vmin=-90,vmax=90)
            ax.set_title(f'{name}\nMAE {metrics[site][name]["mae_deg"]:.1f} deg');ax.set_axis_off()
        fig.colorbar(im,ax=list(axes.flat),label='Axial AoP residual (deg)',shrink=.7)
        role='offset calibration' if site==a.calibration_site else 'fixed-offset check (no site refit)'
        fig.suptitle(f'{site}: {role}; camera basis hypotheses remain provisional')
        fig.savefig(a.output/f'{site}_aop_conventions.png',dpi=160);plt.close(fig)
        if site in DARK_SITES:
            q=d['valid']&d['field']&(abs(d['model'])>.1)
            if q.any():dop[site]=apparent_eta(abs(d['obs'][q]),abs(d['model'][q]),ratios.get(site))
            fig,ax=plt.subplots(figsize=(7,4),layout='constrained');angles=[];measured=[];natural=[]
            for lo in range(10,180,10):
                q=d['valid']&(d['gamma']>=lo)&(d['gamma']<lo+10)
                if q.sum()>=20:
                    angles.append(lo+5);measured.append(float(np.mean(abs(d['obs'][q]))));natural.append(float(np.mean(abs(d['model'][q]))))
            for angle,meas,nat in zip(angles,measured,natural):profile_rows.append(dict(site=site,moon_separation_deg=angle,measured_dop=meas,natural_eta1_dop=nat))
            ax.plot(angles,measured,label='Measured camera DoLP',color='black')
            for eta in [1.,.8,.6]:
                # Missing ALAN deliberately stays unspecified: label natural only.
                addition=ratios.get(site,0.)
                label=f'eta={eta}; '+(f'ALAN/natural={addition}' if site in ratios else 'natural only; ALAN unspecified')
                ax.plot(angles,np.array(natural)*eta/(1+addition),label=label)
            ref=refs[site];ax.set(xlabel='Moon separation (deg)',ylabel='Mean DoLP',ylim=(0,1),title=f'{site}; measured site radiance {float(ref["original_multiyear_monthly_radiance_mean"]):.3f} nW cm^-2 sr^-1')
            ax.legend(fontsize=8);fig.savefig(a.output/f'{site}_dop_atmosphere.png',dpi=170);plt.close(fig)
    with (a.output/'site_radiance_inventory.csv').open('w',newline='') as f:
        keys=list(dict.fromkeys(k for row in inventory for k in row));w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(inventory)
    with (a.output/'aop_metrics.csv').open('w',newline='') as f:
        metric_rows=[dict(site=site,hypothesis=name,offset_calibration_site=a.calibration_site,common_offset_deg=float(np.rad2deg(offsets[name])),**values) for site,rows in metrics.items() for name,values in rows.items()]
        writer=csv.DictWriter(f,fieldnames=list(metric_rows[0]));writer.writeheader();writer.writerows(metric_rows)
    if profile_rows:
        with (a.output/'dark_site_dop_profiles.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(profile_rows[0]));writer.writeheader();writer.writerows(profile_rows)
    sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    result=dict(status='CROSS_SITE_CONVENTION_DIAGNOSTIC_NOT_OPTICAL_CALIBRATION',calibration_site=a.calibration_site,
        common_offsets_deg={k:float(np.rad2deg(v)) for k,v in offsets.items()},aop_metrics=metrics,dark_site_dop=dop,
        reference_sha256=sha(a.reference),anchor_sha256={str(p):sha(p) for p in a.moon_anchors},script_sha256=sha(Path(__file__)),
        capture_alan_ratios_sha256=sha(a.capture_alan_ratios) if a.capture_alan_ratios else None,
        interpretation=['All hypotheses retained; no production convention automatically changed.',
        'Radial spin alternatives diagnose global-versus-local basis mistakes; they are not an instrument calibration.',
        'Cottage estimates offsets; held-out sites use them unchanged. Spatial image parity is not independently established.',
        'Atmosphere factors are prescribed scenarios, not fitted. Unknown ALAN leaves atmospheric eta unidentified.',
        'Skillinge retains its actual radiance even though requested as a dark comparison site.',
        'Historical multiyear upward radiance is not capture-time ground skyglow; never converted directly to an ALAN/natural ratio.',
        'Natural envelope and camera spectral response remain assumptions; no Cottage 0.45 target enters.',
        'Fields are equal-solid-angle equisolid pixels; zenith 30deg radius; audit domain 60deg radius; Moon exclusion 15deg; AoP minimum DoLP 0.1.'])
    (a.output/'convention_and_atmosphere_audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Audit saved to',a.output)


if __name__=='__main__':main()
