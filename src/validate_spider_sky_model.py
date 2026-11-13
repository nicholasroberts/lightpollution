"""Fixed-Moon geometry diagnostics and independent natural-DoLP envelope audit.

Consumes saved measurements only. Never uses adjust_cottage_alan or historical
PRC to choose geometry/amplitude. Unknown camera AoP calibration is explicit.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np


def axial_residual_deg(observed, predicted):
    return np.rad2deg(.5*np.angle(observed*np.conj(predicted)))


def agreement(observed, predicted):
    """Equal solid-angle pixels for an equisolid camera; no DoLP weighting."""
    residual=axial_residual_deg(observed,predicted)
    phase=observed/abs(observed)*np.conj(predicted/abs(predicted))
    offset=.5*np.angle(np.mean(phase))
    centered=axial_residual_deg(observed,predicted*np.exp(2j*offset))
    return dict(pixels=int(len(observed)),raw_axial_mae_deg=float(np.mean(abs(residual))),
        diagnostic_fitted_offset_deg=float(np.rad2deg(offset)),
        offset_removed_axial_mae_deg=float(np.mean(abs(centered))),
        axial_coherence=float(abs(np.mean(phase))))


def envelope(observed,predicted):
    good=np.isfinite(observed)&np.isfinite(predicted)&(predicted>.05)
    ratio=observed[good]/predicted[good]
    if not len(ratio):return dict(pixels=0)
    return dict(pixels=int(len(ratio)),observed_mean=float(np.mean(observed[good])),
        observed_median=float(np.median(observed[good])),model_mean=float(np.mean(predicted[good])),
        ratio_quantiles_10_50_90=np.quantile(ratio,[.1,.5,.9]).tolist(),
        fraction_above_envelope=float(np.mean(ratio>1)),
        mean_absolute_error=float(np.mean(abs(observed[good]-predicted[good]))))


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    from compare_wang2016_measured_skies import compare_site
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',type=Path,required=True)
    p.add_argument('--moon-anchors',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--neutral-separation-deg',type=float,required=True)
    p.add_argument('--reference',type=Path,default=Path(__file__).resolve().parents[1]/'data/reference/prc_spider_18_sites_verified.csv')
    p.add_argument('--sky-masks',type=Path,help='Optional folder of SITE.npy boolean masks in saved Stokes pixel coordinates; True=sky')
    p.add_argument('--low-alan-max-radiance',type=float,default=1.,help='Independent reference VIIRS screen, not proof of ALAN-free sky')
    p.add_argument('--aop-min-dop',type=float,default=.1)
    a=p.parse_args()
    if not 0<a.aop_min_dop<1 or a.low_alan_max_radiance<0:p.error('Invalid threshold')
    if a.output.exists():p.error('Use a new output directory to preserve previous audits')
    a.output.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    anchors=json.loads(a.moon_anchors.read_text())['sites']
    rows=[];results={}
    for ref in csv.DictReader(a.reference.open()):
        site=Path(ref['source_filename']).stem
        row=dict(site=site,site_id=ref['site_id'],status='pending_reviewed_moon_anchor',
            low_alan_candidate=float(ref['original_multiyear_monthly_radiance_mean'])<=a.low_alan_max_radiance)
        rows.append(row)
        if site not in anchors:continue
        directory=a.batch/site/'registration'
        if not (directory/'linear_stokes.npz').exists():row['status']='missing_saved_stokes';continue
        try:
            baseline,arrays=compare_site(directory,None,a.neutral_separation_deg,anchor=anchors[site],alan_relative=0.)
            intensity,obs,rotated,sky,field,center,r90=arrays
            # Undo legacy per-site AoP offset: raw residuals remain available.
            model=rotated*np.exp(-2j*np.deg2rad(baseline['global_aop_offset_deg']))
            y,x=np.indices(obs.shape);r=np.hypot(x-center[0],y-center[1])/r90
            theta=2*np.arcsin(np.clip(r/np.sqrt(2),0,1))
            domain=sky&(theta<=np.deg2rad(60))
            maskpath=a.sky_masks/f'{site}.npy' if a.sky_masks else None
            if maskpath and maskpath.exists():
                mask=np.load(maskpath,allow_pickle=False)
                if mask.dtype!=bool or mask.shape!=obs.shape:raise ValueError('Sky mask must be boolean and match saved Stokes shape')
                domain &= mask
                maskstatus='external_sky_mask_requires_investigator_review'
            else:maskstatus='NO_TERRAIN_MASK_provisional_even_inside_zenith_cap'
            physical=np.isfinite(obs)&(intensity>0)&(abs(obs)<=1)&np.isfinite(model)
            valid=domain&physical
            ao=valid&(abs(obs)>=a.aop_min_dop)&(abs(model)>=a.aop_min_dop)
            if ao.sum()<10:raise ValueError('Fewer than 10 usable AoP pixels')
            stats=dict(site=site,low_alan_candidate=row['low_alan_candidate'],mask_status=maskstatus,
                source=baseline['ephemeris_anchor'],geometry=agreement(obs[ao],model[ao]),
                candidate_elevation_minus_ephemeris_deg=float(90-np.rad2deg(2*np.arcsin(np.linalg.norm(np.asarray(anchors[site]['candidate_pixel_xy'])-center)/(r90*np.sqrt(2))))-anchors[site]['moon_elevation_deg']),
                noncomparable_pixels_in_domain=int((domain&~physical).sum()),
                stokes_sha256=baseline['stokes_sha256'],rotation_report_sha256=baseline['rotation_report_sha256'],
                sky_mask_sha256=sha(maskpath) if maskpath and maskpath.exists() else None,
                natural_envelopes={})
            for region,selection in [('zenith_30deg_radius',valid&field),('extended_60deg_radius',valid)]:
                stats['natural_envelopes'][region]={str(eta):envelope(abs(obs[selection]),eta*abs(model[selection])) for eta in [1.,.8]}
            # Radial scattering-angle profiles compare the predicted DoLP band
            # without maximizing over trees or fitting a field centre.
            anchor=anchors[site];v=np.asarray(anchor['candidate_pixel_xy'])-center
            az=np.arctan2(v[1],v[0]);zs=np.deg2rad(90-anchor['moon_elevation_deg'])
            phi=np.arctan2(y-center[1],x-center[0])
            gamma=np.rad2deg(np.arccos(np.clip(np.cos(theta)*np.cos(zs)+np.sin(theta)*np.sin(zs)*np.cos(phi-az),-1,1)))
            profile=[]
            for lo in range(0,180,10):
                sel=valid&(gamma>=lo)&(gamma<lo+10)
                if sel.sum()>=20:profile.append(dict(scattering_angle_deg=lo+5,pixels=int(sel.sum()),measured=float(np.mean(abs(obs[sel]))),natural=float(np.mean(abs(model[sel])))))
            stats['scattering_profiles']=profile
            fig,axes=plt.subplots(1,4,figsize=(16,4),layout='constrained')
            for ax,data,title,limits,cmap in zip(axes[:3],[abs(obs),abs(model),axial_residual_deg(obs,model)],['Measured DoLP','Natural envelope: eta=1, ALAN=0','Raw AoP residual (uncalibrated)'],[(0,1),(0,1),(-90,90)],['turbo','turbo','twilight']):
                im=ax.imshow(np.where(domain,data,np.nan),vmin=limits[0],vmax=limits[1],cmap=cmap)
                ax.add_patch(plt.Circle(center,r90*np.sin(np.deg2rad(15))/np.sin(np.deg2rad(45)),fill=False,color='white'))
                ax.set_title(title);ax.set_axis_off();fig.colorbar(im,ax=ax,shrink=.65)
            for key,label in [('measured','Measured'),('natural','Natural eta=1')]:axes[3].plot([r['scattering_angle_deg'] for r in profile],[r[key] for r in profile],label=label)
            axes[3].set(xlabel='Moon separation (deg)',ylabel='Mean DoLP',ylim=(0,1));axes[3].legend()
            fig.suptitle(f'{site}: fixed ephemeris + reviewed image bearing; {maskstatus}')
            fig.savefig(a.output/f'{site}_validation.png',dpi=160);plt.close(fig)
            row['status']='provisional_geometry_and_envelope_diagnostics';results[site]=stats
            print(f'{site}: diagnostics saved',flush=True)
        except (ValueError,FileNotFoundError,KeyError) as exc:
            row['status']='blocked';row['reason']=str(exc)
    with (a.output/'site_inventory.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['site','site_id','status','low_alan_candidate','reason']);writer.writeheader();writer.writerows(rows)
    payload=dict(status='diagnostic_not_absolute_polarimetric_validation',sites=results,
        assumptions=['Level equisolid camera; saved nominal analyzer/image-plane basis.',
        'Raw AoP residual uses existing image-plane convention; no full ray-dependent Mueller calibration.',
        'Offset-removed AoP is descriptive shape agreement, not held-out validation.',
        'No fitted Moon position, DoLP amplitude, ALAN contribution or Cottage 0.45 target.',
        'VIIRS low-radiance screen selects candidate low-ALAN sites independently of measured DoLP.',
        'Ratios cannot separate atmosphere, ALAN, calibration and model error.',
        'Measured historical camera channel and green510 animal weighting are not identical spectral measurements.'],
        neutral_separation_deg=a.neutral_separation_deg,aop_min_dop=a.aop_min_dop,
        low_alan_max_radiance=a.low_alan_max_radiance,anchor_sha256=sha(a.moon_anchors),reference_sha256=sha(a.reference),script_sha256=sha(Path(__file__)))
    (a.output/'validation.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    print(f'Wrote {len(results)} site diagnostics; see site_inventory.csv for pending sites.',flush=True)


if __name__=='__main__':main()
