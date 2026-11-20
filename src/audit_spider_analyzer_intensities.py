"""Forward analyzer consistency audit. No fit to astronomical polarization."""
import argparse,json
from itertools import combinations
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from reconstruct_spider_eye_contrast import digest,tiff_channel
from lock_spider_site_registration import warp
from audit_spider_lunar_aop import block_mean


def design(angles):
    a=np.deg2rad(np.asarray(angles))
    return .5*np.column_stack([np.ones(len(a)),np.cos(2*a),np.sin(2*a)])


def solve(data,angles):
    """Unweighted linear Stokes estimate; no gains/offsets fitted."""
    A=design(angles)
    if np.linalg.matrix_rank(A)!=3:raise ValueError('Analyzer design is singular')
    return np.einsum('ij,j...->i...',np.linalg.pinv(A),data)


def forward(stokes,angles):
    return np.einsum('ij,j...->i...',design(angles),stokes)


def leave_one_out(data,angles):
    predictions=[]
    for j in range(4):
        keep=np.arange(4)!=j
        predictions.append(forward(solve(data[keep],np.asarray(angles)[keep]),[angles[j]])[0])
    return np.asarray(predictions)


def stats(values,mask):
    v=np.asarray(values)[...,mask];v=v[np.isfinite(v)]
    if not len(v):return dict(count=0,median=None,median_absolute=None,p95_absolute=None)
    return dict(count=int(len(v)),median=float(np.median(v)),median_absolute=float(np.median(abs(v))),p95_absolute=float(np.percentile(abs(v),95)))


def audit(data,angles,mask):
    stokes=solve(data,angles);prediction=forward(stokes,angles)
    scale=np.mean(data,axis=0)
    with np.errstate(divide='ignore',invalid='ignore'):
        residual=(data-prediction)/scale
        loo=(data-leave_one_out(data,angles))/scale
        dolp=np.hypot(stokes[1],stokes[2])/stokes[0]
    summary=dict(angles_deg=list(map(float,angles)),design_condition_number=float(np.linalg.cond(design(angles))),
        fitted_residual_fraction=stats(residual,mask),leave_one_out_fraction=stats(loo,mask),
        per_analyzer_leave_one_out=[stats(v,mask) for v in loo],
        dolp=stats(dolp,mask),fraction_dolp_above_one=float(np.mean(dolp[mask]>1)))
    return stokes,loo,summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-lock',type=Path,required=True)
    p.add_argument('--site',default='museum')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--bin-size',type=int,default=8)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    if a.bin_size<1:p.error('Bin size must be positive')
    base=a.reference_lock.parent
    checks=json.loads((base/'reference_sha256.json').read_text())
    if checks.get(a.reference_lock.name)!=digest(a.reference_lock):p.error('Registration lock checksum mismatch')
    site=json.loads(a.reference_lock.read_text())['sites'][a.site]
    if site.get('radiometry_status') not in ('black_subtracted_linear_camera_blue','black_subtracted_linear_camera_green'):p.error('Requires linear RAW lock')
    camera_channel=site.get('camera_channel','blue')
    if site.get('status')!='numerically_locked_awaiting_visual_review':p.error('Registration numerical checks must pass')
    groups=site['analyzer_groups']
    sequence_sign=site.get('analyzer_stokes_convention',{}).get('nominal_sequence_sign',1)
    if sequence_sign not in (-1,1):p.error('Invalid signed nominal analyzer convention')
    records={int(Path(r['path']).stem):r for r in site['raw_inputs']}
    transforms={int(Path(r['filename']).stem):r for r in site['transforms']}
    needed=sorted(records)
    data={};input_checks=[]
    for n in needed:
        record=records[n];path=Path(record['path'])
        if digest(path)!=record['sha256']:p.error('Exposure checksum changed: '+str(path))
        image,_=tiff_channel(path,None);t=transforms[n]
        if not t['quality']['numerical_pass']:p.error('Exposure registration failed')
        mapping=t.get('fisheye_mapping')
        if mapping is None:mapping=np.asarray(t['source_to_reference_xy'])
        image=warp(image,mapping) if n!=1 else image
        data[n]=block_mean(image,a.bin_size)
        input_checks.append(dict(exposure=n,path=str(path),sha256=record['sha256']))
    shape=data[needed[0]].shape;y,x=np.indices(shape)
    x=x*a.bin_size+(a.bin_size-1)/2;y=y*a.bin_size+(a.bin_size-1)/2
    c=np.asarray(site['reference_center_xy_px']);r=site['rim_radius_px'];eye=site['eye_radius_px']
    mask=(np.hypot(x-c[0],y-c[1])<eye)
    for v in data.values():mask&=np.isfinite(v)&(v>0)
    if mask.sum()<20:p.error('Insufficient common central samples')
    # One common mask and one common scale; no per-frame gain normalization.
    baseline=float(np.mean(data[needed[0]][mask]))
    sequence=[];within={}
    for n in needed:
        values=data[n][mask];group=(n-1)//3 if max(needed)>4 else n-1
        sequence.append(dict(exposure=n,nominal_angle_deg=group*45,
            datetime_original=records[n].get('datetime_original'),
            mean_linear_camera_channel=float(np.mean(values)),median_linear_camera_channel=float(np.median(values)),
            mean_relative_to_first=float(np.mean(values)/baseline)))
    for group in sorted({(n-1)//3 if max(needed)>4 else n-1 for n in needed}):
        ids=[n for n in needed if ((n-1)//3 if max(needed)>4 else n-1)==group];reference=data[ids[0]]
        within[str(group*45)]=dict(exposures=ids,reference_exposure=ids[0],
            repeat_difference_fraction=[dict(exposure=n,**stats((data[n]-reference)/reference,mask)) for n in ids])
    fig,axes=plt.subplots(2,1,figsize=(12,8),layout='constrained')
    for group,color in enumerate(['tab:blue','tab:orange','tab:green','tab:red']):
        rows=[r for r in sequence if r['nominal_angle_deg']==group*45]
        if not rows:continue
        ids=[r['exposure'] for r in rows]
        axes[0].plot(ids,[r['mean_relative_to_first'] for r in rows],'o-',color=color,label=str(group*45)+' deg')
        values=within[str(group*45)]['repeat_difference_fraction']
        axes[1].plot(ids,[100*r['median'] for r in values],'o-',color=color,label=str(group*45)+' deg')
    axes[0].set_ylabel('Mean intensity / first exposure mean')
    axes[1].set_ylabel('Median pixel change within group (%)')
    for ax in axes:
        ax.set_xticks(range(1,max(needed)+1));ax.set_xlabel('Exposure number in acquisition order');ax.grid(alpha=.2);ax.legend(ncol=4)
    axes[0].set_title('Common central field; common display scale; brightness differences include real polarization')
    axes[1].set_title('Each group relative to its first available exposure; no correction applied')
    fig.suptitle(a.site+': acquisition intensity audit; excluded frames remain absent')
    a.output.mkdir(parents=True)
    fig.savefig(a.output/'exposure_intensity_sequence.png',dpi=150);plt.close(fig)
    results={};maps={};closure=[]
    for k,g in enumerate(groups):
        stack=np.asarray([data[n] for n in g]);scale=stack.mean(0)
        closure.append((stack[0]+stack[2]-stack[1]-stack[3])/(2*scale))
        nominal=[sequence_sign*v for v in [0,45,90,135]]
        rolls=[transforms[n]['image_rotation_clockwise_deg']-transforms[g[0]]['image_rotation_clockwise_deg'] for n in g]
        results[str(g)]={};maps[str(g)]={}
        for label,angles in [('nominal',nominal),('landmark_roll',rolls)]:
            stokes,loo,summary=audit(stack,angles,mask)
            maps[str(g)][label]=(stokes,loo)
            results[str(g)][label]=summary
        results[str(g)]['nominal_pair_sum_closure_fraction']=stats(closure[-1],mask)
    comparison={}
    for ga,gb in combinations(groups,2):
        pair_key=str(ga)+' vs '+str(gb)
        comparison[pair_key]={}
        for label in ['nominal','landmark_roll']:
            s,t=[maps[str(g)][label][0] for g in [ga,gb]]
            with np.errstate(divide='ignore',invalid='ignore'):
                intensity=(s[0]-t[0])/((s[0]+t[0])/2)
                dq=s[1]/s[0]-t[1]/t[0];du=s[2]/s[0]-t[2]/t[0]
                da=.5*np.angle((s[1]+1j*s[2])*np.conj(t[1]+1j*t[2]))*180/np.pi
            comparison[pair_key][label]=dict(intensity_fraction=stats(intensity,mask),normalized_Q_difference=stats(dq,mask),normalized_U_difference=stats(du,mask),aop_difference_deg=stats(da,mask&(np.hypot(s[1],s[2])/s[0]>.05)&(np.hypot(t[1],t[2])/t[0]>.05)))
    fig,axes=plt.subplots(len(groups),6,figsize=(20,4*len(groups)),squeeze=False,layout='constrained')
    for row,g in enumerate(groups):
        stack=np.asarray([data[n] for n in g])
        arrays=[closure[row]]
        titles=['Nominal pair-sum closure']
        for label in ['nominal','landmark_roll']:
            s,loo=maps[str(g)][label]
            arrays.extend([np.sqrt(np.mean(loo**2,axis=0)),np.hypot(s[1],s[2])/s[0]])
            titles.extend([label+': omitted-frame RMS',label+': DoLP'])
        sn=maps[str(g)]['nominal'][0];sr=maps[str(g)]['landmark_roll'][0]
        arrays.append(.5*np.angle((sr[1]+1j*sr[2])*np.conj(sn[1]+1j*sn[2]))*180/np.pi);titles.append('Roll-minus-nominal AoP (deg)')
        for j,(v,title) in enumerate(zip(arrays,titles)):
            ax=axes[row,j];kwargs=dict(cmap='RdBu_r',vmin=-.1,vmax=.1) if j==0 else dict(cmap='viridis',vmin=0,vmax=.2)
            if j in [2,4]:kwargs=dict(cmap='viridis',vmin=0,vmax=1)
            if j==5:kwargs=dict(cmap='RdBu_r',vmin=-10,vmax=10)
            im=ax.imshow(np.where(mask,v,np.nan),**kwargs);ax.set_title(str(g)+'\n'+title,fontsize=9);ax.axis('off');fig.colorbar(im,ax=ax,shrink=.6)
    fig.suptitle(a.site+': analyzer consistency, central field; no sky-pattern, gain or offset fit\nNominal and landmark-roll approximations; fixed analyzer zero unknown. Display scales clipped; JSON statistics unclipped.')
    fig.savefig(a.output/'analyzer_consistency.png',dpi=150);plt.close(fig)
    report=dict(site=a.site,camera_channel=camera_channel,nominal_analyzer_sequence_sign=sequence_sign,status='diagnostic_only_no_correction_applied',reference_lock_sha256=digest(a.reference_lock),script_sha256=digest(Path(__file__)),inputs=input_checks,
        central_sample_count=int(mask.sum()),bin_size=a.bin_size,exposure_sequence=sequence,within_analyzer_groups=within,quartets=results,between_quartets=comparison,
        limitations=['Per-pixel Stokes is estimated from exposures; no astronomical pattern is fitted.',
        'Four samples and three Stokes parameters leave only one residual degree of freedom. Omitted-frame errors are correlated, not independent validation.',
        'Landmark roll is an approximate analyzer angle. Tilt and off-axis rear-analyzer transport require a full optical model.',
        'Unknown constant polarizer axis, transmission, temporal sky changes, vignetting and saturation are not calibrated here.',
        'A sign reversal of all analyzer angles conjugates Q+iU and is intensity-degenerate; intensity closure cannot establish handedness.',
        'Nominal pair sums need not close for non-orthogonal actual analyzer angles.',
        'No theoretical sky rotation, gain, offset, exposure correction or new production Stokes file is applied.'])
    (a.output/'analyzer_consistency.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('Saved',a.output/'analyzer_consistency.png','and analyzer_consistency.json')

if __name__=='__main__':main()
