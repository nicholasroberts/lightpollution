"""Separate spatial intensity weighting from fractional polarization in saved fits."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from compare_spider_mixed_prc import digest,integrated_prc,run as compare


def weighted_pattern(weights,fraction,mask,parallel=.4,perpendicular=.04):
    """Counterfactual Stokes: one intensity map times another fractional field."""
    w=np.asarray(weights)[mask];z=np.asarray(fraction)[mask]
    if not len(w) or not np.isfinite(w).all() or not np.isfinite(z).all() or np.any(w<=0):
        raise ValueError('Invalid common support or weighting')
    d=np.abs(z)
    if np.any(d>1+1e-10):raise ValueError('Fractional polarization exceeds unity')
    mean_d=float(np.sum(w*d)/np.sum(w))
    eff=float(abs(np.sum(w*z))/np.sum(w))
    coherence=float(abs(np.sum(w*z))/np.sum(w*d)) if np.sum(w*d)>0 else None
    result=integrated_prc(w,(w*z).real,(w*z).imag,np.ones(w.shape,bool),parallel,perpendicular)
    result.update(weighted_mean_local_dolp=mean_d,aop_cancellation_factor=coherence,
        integrated_dolp_factorization_error=eff-(mean_d*coherence if coherence is not None else 0))
    return result


def analyze(stokes,mask,parallel=.4,perpendicular=.04):
    fractions={key:(q+1j*u)/i for key,(i,q,u) in stokes.items()}
    weights={key:v[0] for key,v in stokes.items()};weights['uniform']=np.ones(mask.shape)
    matrix={wk:{pk:weighted_pattern(w,z,mask,parallel,perpendicular) for pk,z in fractions.items()}
        for wk,w in weights.items()}
    reference=stokes['measured'][0][mask];reference=reference/reference.mean();intensities={}
    for key,(i,_,_) in stokes.items():
        v=i[mask];v=v/v.mean();delta=v-reference
        correlation=(float(np.corrcoef(v,reference)[0,1]) if np.std(v)>1e-12 and np.std(reference)>1e-12 else None)
        intensities[key]=dict(normalized_rmse=float(np.sqrt(np.mean(delta**2))),
            normalized_median_absolute_difference=float(np.median(abs(delta))),correlation=correlation,
            normalized_p05=float(np.percentile(v,5)),normalized_p95=float(np.percentile(v,95)))
    return matrix,intensities


def run(a):
    if a.output.exists():raise ValueError('Use a fresh output directory')
    # Reuse the existing full provenance, sign, mask and saved-model checks.
    baseline=compare(SimpleNamespace(reference_lock=a.reference_lock,fit=a.fit,output=a.output/'baseline',
        bin_size=a.bin_size,parallel=.4,perpendicular=.04))
    from audit_spider_lunar_aop import block_mean,rays,fixed_moon_in_camera
    from fit_spider_twilight_moon import component,mixture,moon_plane_reflection
    fit=json.loads((a.fit/'twilight_moon_fit.json').read_text());site=fit['site']
    lock=json.loads(a.reference_lock.read_text())['sites'][site]
    with np.load(a.reference_lock.parent/site/'registration/linear_stokes.npz',allow_pickle=False) as z:
        I,Q,U=[block_mean(z[k],a.bin_size) for k in ['I','Q','U']]
    yy,xx=np.indices(I.shape);xy=np.stack([xx*a.bin_size+(a.bin_size-1)/2,yy*a.bin_size+(a.bin_size-1)/2],axis=-1)
    c=np.asarray(lock['reference_center_xy_px']);view=rays(xy,c,lock['rim_radius_px'])
    m=fixed_moon_in_camera(lock['astronomical_anchor'],lock['camera_to_north_xy'])
    solar=fit['solar'];s=fixed_moon_in_camera(dict(moon_elevation_deg=solar['sun_elevation_deg'],
        moon_compass_azimuth_deg=solar['sun_compass_azimuth_deg']),lock['camera_to_north_xy'])
    if fit['handedness_hypothesis']=='reversed':
        reflection=moon_plane_reflection(c,lock['moon_pixel_xy'])
        m=np.r_[reflection@m[:2],m[2]];s=np.r_[reflection@s[:2],s[2]]
    m=component(view,m,fit['fixed_pmax']);s=component(view,s,fit['fixed_pmax'])
    models=dict(moon=(m[0]+fit['models']['moon']['parameters'][0],m[1]),
        sun=(s[0]+fit['models']['sun']['parameters'][0],s[1]),mixed=mixture(m,s,*fit['models']['mixed']['parameters']))
    stokes=dict(measured=(I,Q,U),**{k:(i,z.real,z.imag) for k,(i,z) in models.items()})
    with np.load(a.output/'baseline/prc_support.npz',allow_pickle=False) as z:regions={k:np.array(z[k]) for k in z.files}
    results={}
    for region,mask in regions.items():
        matrix,stats=analyze(stokes,mask)
        for key in stokes:
            if not np.isclose(matrix[key][key]['maximum_absolute_prc'],baseline['results'][region][key]['maximum_absolute_prc'],rtol=1e-12,atol=1e-12):
                raise ValueError('Native weighting does not reproduce baseline')
        results[region]=dict(weight_by_pattern=matrix,normalized_intensity_comparison=stats)
    report=dict(site=site,camera_channel=fit['camera_channel'],results=results,
        common_circle_coverage_fraction=baseline['common_circle_coverage_fraction'],
        matrix_definition='Rows supply I weights; columns supply (Q+iU)/I. Diagonal reproduces original PRC. Other entries are counterfactual diagnostics.',
        factorization='integrated DoLP = intensity-weighted mean local DoLP * AoP cancellation factor',
        cancellation_definition='abs(sum(I * DoLP * exp(2i AoP))) / sum(I * DoLP); 1=no cancellation',
        source_lock_sha256=baseline['source_lock_sha256'],source_stokes_sha256=baseline['source_stokes_sha256'],
        source_fit_sha256=baseline['source_fit_sha256'],source_fit_arrays_sha256=baseline['source_fit_arrays_sha256'],
        script_sha256=digest(__file__),baseline_script_sha256=baseline['script_sha256'],
        fit_or_zero_calibration_applied=False,limitations=baseline['limitations']+[
        'Swapping normalized intensity weights is a diagnostic, not a corrected physical Sun/Moon model.',
        'Intensity maps are normalized on the common support; absolute source radiometry is not tested.',
        'DoLP and AoP contributions interact; the factorization does not uniquely assign additive causal error shares.'])
    (a.output/'prc_weighting_audit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    keys=list(stokes);weightkeys=keys+['uniform'];mask=regions['common_central_field']
    matrix=results['common_central_field']['weight_by_pattern']
    fig,axes=plt.subplots(2,4,figsize=(14,7),layout='constrained')
    for ax,key in zip(axes[0],keys):
        i=stokes[key][0];v=np.where(mask,i/np.mean(i[mask]),np.nan)
        im=ax.imshow(v,cmap='viridis',vmin=.7,vmax=1.3);ax.set_title(key+': I / field mean');ax.axis('off')
        ys,xs=np.where(mask);ax.set_xlim(xs.min()-2,xs.max()+2);ax.set_ylim(ys.max()+2,ys.min()-2)
        fig.colorbar(im,ax=ax,shrink=.7)
    im=axes[1,0].imshow([[matrix[w][p]['maximum_absolute_prc'] for p in keys] for w in weightkeys],cmap='viridis')
    axes[1,0].set_xticks(range(4));axes[1,0].set_xticklabels(keys,rotation=30)
    axes[1,0].set_yticks(range(5));axes[1,0].set_yticklabels(weightkeys)
    axes[1,0].set_xlabel('Fractional polarization pattern');axes[1,0].set_ylabel('Intensity weights')
    axes[1,0].set_title('Maximum |PRC|');fig.colorbar(im,ax=axes[1,0],shrink=.7)
    for ax,metric,title in zip(axes[1,1:],['maximum_absolute_prc','weighted_mean_local_dolp','aop_cancellation_factor'],
            ['Maximum |PRC|','Weighted local DoLP','AoP cancellation factor']):
        for w in ['measured','uniform']:
            ax.plot(range(4),[matrix[w][p][metric] for p in keys],'o-',label=w+' I weights')
        ax.set_xticks(range(4));ax.set_xticklabels(keys,rotation=30);ax.set_title(title);ax.legend(fontsize=8)
    fig.suptitle(site+' — common central field; fixed fits; no analyzer-zero adjustment\nIntensity displays clipped to 0.7–1.3; numerical statistics unclipped')
    fig.savefig(a.output/'prc_weighting_audit.png',dpi=160);plt.close(fig)
    print('Saved',a.output/'prc_weighting_audit.json')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-lock',type=Path,required=True)
    p.add_argument('--fit',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--bin-size',type=int,default=8)
    a=p.parse_args()
    if a.bin_size<1:p.error('Bin size must be positive')
    try:run(a)
    except (OSError,ValueError,KeyError) as exc:p.error(str(exc))


if __name__=='__main__':main()
