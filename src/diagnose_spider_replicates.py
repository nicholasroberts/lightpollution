"""Diagnose individual sky replicates using saved transforms, without registration.

Read original TIFFs, verify hashes and apply exactly the saved transforms.
Preserve nominal analyzer order. Produce per-replicate and combined summaries,
DoP/AoP maps and eye-axis angles, with fixed 30deg-radius integration.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from reconstruct_spider_eye_contrast import (digest,tiff_channel,circle_mask,
    combine_replicates,eye_result,GROUPS,reference_comparison)
from analyse_spider_sky_rotation import warp_to_reference
from audit_spider_polarization import summarize,original_aop_cmap


def replicate_results(quartets,mask):
    """Same calculation and same mask for each replicate; no fitting to reference."""
    results={}
    for method in ['linear_stokes','legacy_dop_aop']:
        entries=[]
        for label,selection in [(f'replicate_{n+1}',[q]) for n,q in enumerate(quartets)]+[('combined',quartets)]:
            i,q,u,qa=combine_replicates(iter(selection),method)
            stats,_,_=summarize(i,q,u,mask)
            for mode in ['stokes_intensity','normalized_direction']:
                result=eye_result(i,q,u,mask,response_mode=mode,invalid_policy='exclude')
                stats['fixed_eye_results'][mode]['max_axis_deg_mod_180']=result['max_axis_deg_mod_180']
                stats['fixed_eye_results'][mode]['sum_I']=result['sum_I']
                stats['fixed_eye_results'][mode]['sum_Q']=result['sum_Q']
                stats['fixed_eye_results'][mode]['sum_U']=result['sum_U']
            stats['label']=label
            stats['pair_closure']=[dict(median_abs=float(np.nanmedian(abs(r['closure'][mask]))),
                p95_abs=float(np.nanpercentile(abs(r['closure'][mask]),95))) for r in qa]
            entries.append(stats)
        results[method]=entries
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();source=a.run.resolve();out=a.output.resolve()
    if out.exists():p.error('Use a new output directory')
    if out==source or source in out.parents or out in source.parents:p.error('Output must be separate from saved run')
    report_path=source/'rotation_analysis.json'
    report=json.loads(report_path.read_text())
    inputs={Path(r['path']).name:r for r in report['inputs']}
    numbers={int(Path(n).stem) for n in inputs}
    groups=GROUPS if numbers==set(range(1,13)) else [[1,2,3,4]] if numbers=={1,2,3,4} else None
    if groups is None:p.error('Expected complete 12-image set or one four-image quartet')
    transforms={r['filename']:r for r in report['rotations']}
    center=report['reference_center_xy_px']
    radius=report.get('eye_radius_px',report['radius_60_px']*np.sin(np.deg2rad(15))/np.sin(np.deg2rad(30)))
    if report.get('eye_angular_radius_deg',30)!=30:p.error('Saved report has incompatible eye radius')
    aligned={};input_qa=[]
    for number in sorted(numbers):
        name=f'{number}.tif';record=inputs[name];transform=transforms[name]
        print(f'Applying saved transform: {name}',flush=True)
        path=Path(record['path'])
        if path.resolve()==out or path.resolve() in out.parents:p.error('Output overlaps raw file')
        image,metadata=tiff_channel(path,record['channel'])
        if metadata['sha256']!=record['sha256']:raise ValueError(f'Source hash mismatch: {name}')
        refinement=transform['refinement']
        aligned[number]=warp_to_reference(image,refinement['refined_source_center_xy_px'],center,
                                         transform['scene_roll_clockwise_deg'])
        input_qa.append(metadata)
    quartets=[[aligned[n] for n in g] for g in groups]
    mask=circle_mask(quartets[0][0].shape,center,radius)
    results=replicate_results(quartets,mask)
    reference=reference_comparison(report['site_id'],results['linear_stokes'][-1]['fixed_eye_results']['stokes_intensity']['prc'])
    expected=reference['reference_rounded_prc']
    for entries in results.values():
        for entry in entries:
            for result in entry['fixed_eye_results'].values():
                result['difference_from_original']=result['prc']-expected
    out.mkdir(parents=True)
    for method in results:
        fig,axes=plt.subplots(len(groups)+1,3,figsize=(12,3.6*(len(groups)+1)),squeeze=False)
        selections=[[q] for q in quartets]+[quartets]
        for row,selection in enumerate(selections):
            i,q,u,_=combine_replicates(iter(selection),method)
            with np.errstate(divide='ignore',invalid='ignore'):
                dop=np.hypot(q,u)/i;angle=np.rad2deg(.5*np.arctan2(u,q))
            vmax=float(np.nanpercentile(i[mask],99))
            for col,(arr,label,cm,lo,hi) in enumerate([
                (i,'Intensity','gray',0,vmax),(dop,'DoLP','turbo',0,1),
                (angle,'AoP (degrees)',original_aop_cmap(),-90,90)]):
                ax=axes[row,col]
                if isinstance(cm,str):cm=plt.get_cmap(cm).copy()
                cm.set_bad('lightgray')
                artist=ax.imshow(arr.T,cmap=cm,vmin=lo,vmax=hi)
                ax.add_patch(Circle((center[1],center[0]),radius,fill=False,color='black' if col else 'yellow',lw=.8))
                # Crop display to near the eye field for readable detail; mask not altered.
                ax.set_xlim(center[1]-1.3*radius,center[1]+1.3*radius)
                ax.set_ylim(center[0]+1.3*radius,center[0]-1.3*radius)
                result=results[method][row]['fixed_eye_results']['stokes_intensity']
                ax.set_title(f'{results[method][row]["label"]}: {label}\nPRC {result["prc"]:.4f}')
                ax.set_aspect('equal');ax.axis('off');fig.colorbar(artist,ax=ax,shrink=.7)
        fig.suptitle(f'{report["site_id"]}: {method}; saved alignment reused\n60deg full field; 30deg radius; no new registration or parameter fitting')
        fig.tight_layout(rect=[0,0,1,.95]);fig.savefig(out/f'{method}_replicates.png',dpi=160);plt.close(fig)
    payload=dict(status='diagnostic_not_historical_validation',site_id=report['site_id'],reference=reference,
        source_report_sha256=digest(report_path),script_sha256=digest(__file__),input_files=input_qa,
        full_field_deg=60,angular_radius_deg=30,center_xy_px=center,radius_px=radius,
        nominal_groups=groups,methods=results,raw_images_modified=False,
        limitations=['Saved residual alignment errors remain unchanged',
        'Individual eye-axis angles are relative to common image coordinates, not geographic bearings',
        'Legacy arithmetic AoP averaging retained only as explicit comparison',
        'No tuning or method selection by agreement with reference'])
    (out/'replicate_summary.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    print(f'Wrote replicate_summary.json and two replicate-map PNGs in {out}',flush=True)

if __name__=='__main__':main()
