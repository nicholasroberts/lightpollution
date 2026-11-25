"""Figure 1B: four lunar-only maps forming one subfigure; no recalculation."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
try:
    from cmcrameri import cm as scientific_cmaps
except ImportError as exc:
    raise SystemExit('Install figure palettes: python -m pip install cmcrameri') from exc
from reconstruct_spider_eye_contrast import digest


def align_analyzer_basis(model_q, model_u, measured_q, measured_u, mask, withheld):
    """One constant spin-2 phase; no spatial rotation or change to DoLP."""
    model=model_q+1j*model_u;measured=measured_q+1j*measured_u
    good=mask & np.isfinite(model) & np.isfinite(measured) & (abs(model)>1e-10) & (abs(measured)>1e-10)
    train=good & ~withheld;test=good & withheld
    if min(train.sum(),test.sum())<20:raise ValueError('Insufficient analyzer alignment support')
    mean=np.mean(measured[train]/abs(measured[train])*np.conj(model[train]/abs(model[train])))
    if abs(mean)<1e-6:raise ValueError('Undefined analyzer offset')
    offset=.5*np.angle(mean);aligned=model*np.exp(2j*offset)
    metrics={}
    for name,selection in [('training',train),('withheld',test)]:
        error=np.degrees(.5*np.angle(measured[selection]*np.conj(aligned[selection])))
        metrics[name]=dict(pixel_count=int(selection.sum()),rms_deg=float(np.sqrt(np.mean(error**2))))
    return np.degrees(.5*np.angle(aligned)),dict(model_to_measured_offset_deg=float(np.degrees(offset)),
        training_circular_coherence=float(abs(mean)),metrics=metrics,
        method='Equal-weight mean unit spin-2 phase difference on training pixels; checkerboard 96px blocks, every third partition withheld; central 60deg field, measured DoLP 0.05 to 0.95.',
        calibration_status='Pattern-estimated constant analyzer zero, not independent instrument calibration',
        spatial_rotation_fitted=False,DoLP_changed=False,
        interpretation='Lunar-only conditional basis alignment; museum structural disagreement is retained. Combined-model alignment must be estimated once from its combined pattern and shared by all components.')

def load_site(folder, lock_path=None, north=None):
    rp=folder/'lunar_only.json';ap=folder/'lunar_only_stokes.npz'
    report=json.loads(rp.read_text())
    with np.load(ap,allow_pickle=False) as z:
        d=np.array(z['DoLP']);angle=np.array(z['AoP_deg'])
        if 'central_60deg_field' in z:
            field=np.array(z['central_60deg_field'])
            c=np.array([(d.shape[1]-1)/2,(d.shape[0]-1)/2])
            radius=(d.shape[0]-1)/2-2
            eye=radius*np.sqrt(2)*np.sin(np.pi/12)
        else:
            if lock_path is None:
                raise ValueError('Source lock required to reconstruct complete theoretical field')
            if report['source_lock_sha256']!=digest(lock_path):
                raise ValueError('Source lock changed')
            lock=json.loads(lock_path.read_text())['sites'][report['site']]
            c=np.asarray(lock['reference_center_xy_px'])
            radius=lock['rim_radius_px'];eye=lock['eye_radius_px']
            y,x=np.indices(d.shape);field=np.hypot(x-c[0],y-c[1])<eye
    valid=field&np.isfinite(d)&(d>=0)&(d<=1)
    if valid.sum()<2:
        raise ValueError('Insufficient theoretical field values')
    row=dict(site=report['site'],mean_dolp=float(d[valid].mean()),
        spatial_sd=float(d[valid].std(ddof=1)),included_pixels=int(valid.sum()),
        excluded_pixels=int((field&~valid).sum()),support='Full theoretical zenith-centred 60-degree field; equal pixels; sample SD ddof=1')
    anchor=report.get('astronomical_anchor',report)
    for key in ['local_time','utc','moon_elevation_deg','moon_compass_azimuth_deg']:
        row[key]=anchor.get(key)
    if report.get('neutral_separation_deg')!=70:
        raise ValueError('This figure specification requires shared 70-degree separation')
    row['neutral_separation_deg']=report['neutral_separation_deg']
    row['source_report_sha256']=digest(rp);row['source_arrays_sha256']=digest(ap)
    if north is not None:
        npz=north/'north_stokes.npz';nr=north/'north_rotation.json'
        if digest(npz)!=report['source_north_arrays_sha256'] or digest(nr)!=report['source_north_report_sha256']:
            raise ValueError('Measured north arrays/report differ from model inputs')
        with np.load(npz,allow_pickle=False) as measured, np.load(ap,allow_pickle=False) as model:
            md=measured['DoLP'];support=field & measured['valid'] & np.isfinite(md) & (md>=.05) & (md<=.95)
            yy,xx=np.indices(d.shape);withheld=(xx//96+yy//96)%3==0
            aligned,registration=align_analyzer_basis(model['Q'],model['U'],measured['Q'],measured['U'],support,withheld)
        angle=np.where(np.isfinite(angle),aligned,np.nan)
        row['analyzer_basis_alignment']=registration
        row['alignment_source_north_sha256']=digest(npz)
    return report,d,angle,c,radius,eye,row


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cottage',type=Path,required=True)
    p.add_argument('--museum',type=Path,required=True)
    p.add_argument('--museum-lock',type=Path,required=True)
    p.add_argument('--cottage-lock',type=Path,help='Required for Cottage theory calculated from its RAW celestial lock')
    p.add_argument('--cottage-north',type=Path,required=True)
    p.add_argument('--museum-north',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    panels=[load_site(a.cottage,a.cottage_lock,a.cottage_north),load_site(a.museum,a.museum_lock,a.museum_north)]
    if [v[0]['site'] for v in panels]!=['Cottage','museum']:
        p.error('Expected Cottage then museum reports')
    if panels[0][0]['kernel_sha256']!=panels[1][0]['kernel_sha256']:
        p.error('Berry/Wang kernel versions differ; reproduce both with the same kernel')
    plt.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,2,figsize=(9,9),layout='constrained')
    for ri,(report,d,angle,c,r,eye,row) in enumerate(panels):
        el,az=np.deg2rad([row['moon_elevation_deg'],row['moon_compass_azimuth_deg']])
        mr=r*np.sqrt(1-np.sin(el))
        moon=c+np.array([-np.sin(az),-np.cos(az)])*mr
        for col,(v,cmap,lo,hi,label) in enumerate([(d,plt.get_cmap('viridis'),0,1,'DoLP'),(angle,scientific_cmaps.romaO,-90,90,'AoP')]):
            ax=axes[ri,col];im=ax.imshow(v,cmap=cmap,vmin=lo,vmax=hi,rasterized=True)
            ax.add_patch(plt.Circle(c,eye,fill=False,color='white',lw=.9))
            ax.plot(*moon,'r+',ms=6)
            for tx,ty,t in [(.5,.99,'N'),(.01,.5,'E'),(.5,.01,'S'),(.99,.5,'W')]:
                ax.text(tx,ty,t,transform=ax.transAxes,ha='center',va='center',fontsize=9)
            title=report['site']+' — '+label
            if col==0:title+='\n'+f"{row['mean_dolp']:.3f} ± {row['spatial_sd']:.3f} (spatial SD)"
            ax.set_title(title,fontsize=11)
            ax.set_xlim(c[0]-1.07*r,c[0]+1.07*r);ax.set_ylim(c[1]+1.07*r,c[1]-1.07*r)
            ax.axis('off')
    fig.colorbar(axes[0,0].images[0],ax=axes[:,0].tolist(),shrink=.65,label='DoLP')
    fig.colorbar(axes[0,1].images[0],ax=axes[:,1].tolist(),shrink=.65,label='AoP (degrees modulo 180)')
    fig.suptitle('Lunar-only Berry/Wang baseline\n70° separation; no ALAN or dilution; constant analyzer zero aligned per site',fontsize=13)
    fig.text(.015,.985,'B',ha='left',va='top',fontsize=18,fontweight='bold')
    a.output.mkdir(parents=True)
    for ext in ['png','pdf','svg']:fig.savefig(a.output/('fig01b_two_sites.'+ext),dpi=300)
    plt.close(fig)
    rows=[v[-1] for v in panels]
    with (a.output/'figure_values.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (a.output/'figure_values.json').write_text(json.dumps(dict(sites=rows,script_sha256=digest(__file__),colormaps=dict(DoLP='viridis',AoP='romaO'),dolp_colour_normalization='linear',
        definition='Theoretical only; complete 60-degree full field, not measured-mask intersection. Numerical means may differ slightly from prior museum common-support diagnostics. No animal spectral claim. AoP uses one pattern-estimated constant analyzer zero per site, not an independently calibrated absolute angle.'),indent=2)+'\n')
    lines=['# Figure 1B theoretical field statistics','',
        'Mean of pixel DoLP ± 1 spatial sample SD; zenith-centred 60° full field.',
        'These are theoretical predictions, not experimental measurements.','',
        '| Site | Mean DoLP | Spatial SD | Included | Excluded |',
        '| --- | --- | --- | --- | --- |']
    lines += [f"| {r['site']} | {r['mean_dolp']:.9f} | {r['spatial_sd']:.9f} | {r['included_pixels']} | {r['excluded_pixels']} |" for r in rows]
    (a.output/'labbook_values.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(rows,indent=2))
    print('Saved',a.output)


if __name__=='__main__':main()
