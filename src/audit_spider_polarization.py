"""Compare existing reconstruction arrays; no re-registration or parameter fitting.

Matched scales: DoLP 0..1, AoP -90..90 degrees. Transpose option matches R's
x=image-row/y=image-column display only; scientific masks remain in TIFF axes.
The scalar diagnostics distinguish local polarization from cancellation across
sky directions. A local-max benchmark is explicitly not a pair of fixed eyes.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from reconstruct_spider_eye_contrast import circle_mask,eye_result,digest


# Approximate colour progression sampled from investigator's original PNG legend.
# Its numeric unit/range is unverified; do not infer or modify angles from it.
ORIGINAL_AOP_COLORS=['#ff2500','#ff7900','#e6f001','#5ef901','#02f943',
    '#04fcda','#008dff','#1933ff','#a038ff','#ff37ce','#ff284d']

def original_aop_cmap():
    return LinearSegmentedColormap.from_list('original_aop_legend',ORIGINAL_AOP_COLORS,N=256)


def quantiles(values):
    values=values[np.isfinite(values)]
    if not len(values): return None
    return dict(zip(['min','p05','p25','median','p75','p95','max'],
        map(float,np.quantile(values,[0,.05,.25,.5,.75,.95,1]))))


def summarize(i,q,u,mask):
    with np.errstate(divide='ignore',invalid='ignore'):
        dop=np.hypot(q,u)/i; angle=.5*np.arctan2(u,q)
    finite=mask & np.isfinite(dop) & np.isfinite(angle) & (i>0)
    physical=finite & (dop<=1+1e-10)
    d=dop[physical]; a=angle[physical]; iv=i[physical]
    if not len(d): raise ValueError('No physical FOV samples')
    coherence=abs(np.mean(d*np.exp(2j*a)))/np.mean(d) if np.mean(d)>0 else None
    ratios={}
    for mode in ['stokes_intensity','normalized_direction']:
        result=eye_result(i,q,u,mask,response_mode=mode,invalid_policy='exclude')
        ratios[mode]={k:result[k] for k in ['prc','effective_integrated_DoLP','excluded_pixels','excluded_fraction']}
    # Diagnostic only: replace orientation differences by local alignment, which
    # cannot represent two fixed polarization axes over the whole field.
    absorption_p=10**(-.4); absorption_t=10**(-.04)
    b=1-(absorption_p+absorption_t)/2; gain=(absorption_t-absorption_p)/2
    mean_d=float(np.mean(d)); weighted_d=float(np.sum(iv*d)/np.sum(iv))
    benchmark=lambda degree: float(np.log((b+gain*degree)/(b-gain*degree)))
    summary=dict(fov_pixels=int(mask.sum()),physical_pixels=int(physical.sum()),
        nonfinite_or_nonpositive_pixels=int((mask & ~finite).sum()),
        dop_above_one_pixels=int((finite & (dop>1+1e-10)).sum()),
        intensity_quantiles=quantiles(i[mask]),dop_quantiles=quantiles(dop[finite]),
        physical_mean_dop=mean_d,intensity_weighted_mean_dop=weighted_d,
        orientation_coherence_unweighted=None if coherence is None else float(coherence),
        # AoP is axial: resultant uses twice the angle; no arithmetic mean used.
        fixed_eye_results=ratios,
        hypothetical_locally_aligned_benchmark=dict(
            warning='Not the fixed two-eye model; diagnostic for angular cancellation only',
            normalized_direction_prc=benchmark(mean_d),stokes_intensity_prc=benchmark(weighted_d)))
    return summary,dop,np.rad2deg(angle)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--angular-radius-deg',type=float,default=30.,help='Eye angular radius; default 30 means 60deg full field')
    a=p.parse_args(); source=a.run.resolve(); out=a.output.resolve()
    if out==source or source in out.parents or out in source.parents:
        p.error('Output must be a separate sibling directory')
    if out.exists(): p.error('Output already exists')
    report=json.loads((source/'rotation_analysis.json').read_text())
    if not 0<a.angular_radius_deg<=90: p.error('Angular radius must be in (0,90]')
    center=report['reference_center_xy_px']
    # Existing registration report radius refers to the earlier 60deg RADIAL mask.
    radius=report['radius_60_px']*np.sin(np.deg2rad(a.angular_radius_deg/2))/np.sin(np.deg2rad(30))
    field_label=f'{2*a.angular_radius_deg:g}deg full field ({a.angular_radius_deg:g}deg radius)'
    fig,axes=plt.subplots(2,3,figsize=(13,9))
    results={}
    aop_fig,aop_axes=plt.subplots(1,2,figsize=(14,8))
    for row,method in enumerate(['linear_stokes','legacy_dop_aop']):
        path=source/f'{method}.npz'
        with np.load(path,allow_pickle=False) as z: i,q,u=[z[k] for k in ['I','Q','U']]
        mask=circle_mask(i.shape,center,radius)
        stats,dop,angle=summarize(i,q,u,mask)
        aa=aop_axes[row]
        cmap=original_aop_cmap();cmap.set_bad('lightgray')
        art=aa.imshow(angle.T,cmap=cmap,vmin=-90,vmax=90,aspect='equal')
        aa.contour(mask.T,levels=[.5],colors='black',linewidths=.7)
        aa.set_title(method);aa.set_xlabel('Image row (original R x)');aa.set_ylabel('Image column (original R y)')
        cb=aop_fig.colorbar(art,ax=aa,shrink=.8,ticks=[-90,-45,0,45,90]);cb.set_label('AoP (degrees)')
        stats['source_sha256']=digest(path);results[method]=stats
        # Full sky in column 1, then cropped physical field in both orientations.
        arrays=[dop,np.where(mask,dop,np.nan).T,np.where(mask,angle,np.nan).T]
        titles=[f'{method}: full DoLP', f'DoLP: {field_label}',f'AoP: {field_label}']
        for col,(ax,arr,title) in enumerate(zip(axes[row],arrays,titles)):
            cm=plt.get_cmap('turbo').copy() if col<2 else original_aop_cmap();cm.set_bad('lightgray')
            artist=ax.imshow(arr,cmap=cm,vmin=0 if col<2 else -90,vmax=1 if col<2 else 90)
            ax.set_title(title);ax.set_aspect('equal');ax.axis('off');fig.colorbar(artist,ax=ax,shrink=.7)
    fig.suptitle('MATCHED POLARIZATION DIAGNOSTICS — fixed display scales; arrays unchanged\nTransposed panels match R axis ordering only; absolute AoP zero/sign not established')
    fig.tight_layout(rect=[0,0,1,.93])
    out.mkdir(parents=True)
    fig.savefig(out/'matched_polarization.png',dpi=170);plt.close(fig)
    aop_fig.suptitle(f'AoP comparison — palette approximates original PNG legend; shared -90 to +90 degrees\nBlack circle = {field_label}; original legend units remain unverified')
    aop_fig.tight_layout(rect=[0,0,1,.93]);aop_fig.savefig(out/'matched_aop.png',dpi=180);plt.close(aop_fig)
    result=dict(aop_display=dict(palette='original_PNG_legend_approximation',limits_degrees=[-90,90],
        underlying_angles_modified=False),status='diagnostic_not_historical_validation',
        analysis_sha256=digest(source/'rotation_analysis.json'),script_sha256=digest(__file__),
        center_xy_px=center,angular_radius_deg=a.angular_radius_deg,full_field_deg=2*a.angular_radius_deg,
        radius_px=radius,previous_60deg_radial_radius_px=report['radius_60_px'],methods=results,
        limitations=['Physical-pixel exclusion is reported, not hidden',
        'Original PNG colours are not numeric ground truth',
        'Transpose is display-only; AoP basis is unchanged',
        'Recovered R script converts AoP to degrees; original figure legend units remain unverified'],
        raw_images_modified=False)
    (out/'polarization_summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(f'Wrote {out}/matched_polarization.png and polarization_summary.json')

if __name__=='__main__': main()
