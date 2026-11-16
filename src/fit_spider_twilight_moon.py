"""Exploratory fixed-geometry DoLP fit; ideal solar/lunar scattering proxies."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from audit_spider_lunar_aop import (digest,block_mean,rays,lunar_phase,
    fixed_moon_in_camera,solar_geometry,residual_summary,pixels)
from lock_spider_site_registration import north_rotation


def moon_plane_reflection(center,moon_pixel):
    """Orthogonal reflection about observed Moon/zenith meridian, not a fit."""
    u=np.asarray(moon_pixel,dtype=float)-center
    if np.linalg.norm(u)<1:raise ValueError('Moon too close to zenith to anchor handedness')
    u=u/np.linalg.norm(u)
    return 2*np.outer(u,u)-np.eye(2)


def component(view,source,pmax):
    """Rayleigh-shaped intensity and complex linear polarization, zenith I=1.

    A geometric proxy, not twilight radiative transfer or the production KS91 model.
    """
    mu=np.sum(view*source,axis=-1)
    norm=1+source[2]**2
    intensity=(1+mu**2)/norm
    phase=lunar_phase(view,source)
    polarized=pmax*(1-mu**2)/norm
    # At exact source/antisource the polarized contribution is zero.
    z=np.where(polarized<1e-12,0,polarized*phase)
    return intensity,z


def mixture(moon,sun,ratio,background):
    im,zm=moon;iss,zs=sun
    return im+ratio*iss+background,zm+ratio*zs


def fit_model(observed,moon,sun,train,kind):
    def prediction(q):
        if kind=='moon':return mixture(moon,sun,0,q[0])
        if kind=='sun':return sun[0]+q[0],sun[1]
        return mixture(moon,sun,q[0],q[1])
    def residual(q):
        intensity,z=prediction(q)
        return (abs(z)/intensity-observed)[train]
    starts=[[0],[1],[10]] if kind!='mixed' else [[r,b] for r in [.01,.1,1,10,100] for b in [0,1,10]]
    bounds=([0],[1000]) if kind!='mixed' else ([0,0],[1000,1000])
    fits=[least_squares(residual,q,bounds=bounds,loss='linear',max_nfev=1000) for q in starts]
    best=min(fits,key=lambda r:np.sum(r.fun**2))
    intensity,z=prediction(best.x)
    return dict(parameters=best.x.tolist(),success=bool(best.success),cost=float(best.cost),
        upper_bound_reached=bool(np.any(best.x>999)),jacobian_condition_number=float(np.linalg.cond(best.jac))),intensity,z


def errors(obs,pred,mask):
    d=(pred-obs)[mask]
    return dict(count=int(len(d)),rmse=float(np.sqrt(np.mean(d*d))),median_absolute_error=float(np.median(abs(d))),mean_bias=float(np.mean(d)))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-lock',type=Path,required=True)
    p.add_argument('--site',default='museum')
    p.add_argument('--handedness',choices=['original','reversed'],default='original',help='Explicit camera parity hypothesis; no automatic selection')
    p.add_argument('--reference',type=Path,default=Path(__file__).resolve().parents[1]/'data/reference/prc_spider_18_sites_verified.csv')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--bin-size',type=int,default=8)
    p.add_argument('--pmax',type=float,default=.8,help='Fixed ideal maximum DoLP for each source; not a fitted atmospheric factor')
    p.add_argument('--sky-mask',type=Path,help='Optional original-shape NPY boolean clear-sky mask; otherwise fit central 60deg full field')
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    if a.bin_size<1 or not 0<a.pmax<=1:p.error('Invalid bin size or pmax')
    base=a.reference_lock.parent;npz=base/a.site/'registration/linear_stokes.npz'
    checks=json.loads((base/'reference_sha256.json').read_text())
    for f in [a.reference_lock,npz]:
        if checks.get(str(f.relative_to(base)))!=digest(f):p.error('Locked input checksum mismatch: '+str(f))
    site=json.loads(a.reference_lock.read_text())['sites'][a.site]
    if site.get('radiometry_status')!='black_subtracted_linear_camera_blue':p.error('Requires the linear RAW registration lock')
    if site.get('status')!='numerically_locked_awaiting_visual_review':p.error('Registration numerical checks must pass')
    with np.load(npz) as z:
        original_shape=z['I'].shape
        I,Q,U=[block_mean(z[k],a.bin_size) for k in ['I','Q','U']]
    y,x=np.indices(I.shape);x=x*a.bin_size+(a.bin_size-1)/2;y=y*a.bin_size+(a.bin_size-1)/2
    c=np.asarray(site['reference_center_xy_px']);r=site['rim_radius_px']
    if site.get('moon_pixel_xy') is None:p.error('Locked observed Moon pixel is required')
    observed_moon=np.asarray(site['moon_pixel_xy'])
    north_matrix,north_angle=north_rotation(c,observed_moon,site['astronomical_anchor']['moon_compass_azimuth_deg'])
    if not np.allclose(north_matrix,site['camera_to_north_xy'],atol=1e-6):p.error('Saved north rotation disagrees with locked observed Moon; audit registration lock first')
    view=rays(np.stack([x,y],axis=-1),c,r)
    m=fixed_moon_in_camera(site['astronomical_anchor'],site['camera_to_north_xy'])
    s,solar=solar_geometry(site['astronomical_anchor'],a.site,a.reference,site['camera_to_north_xy'])
    if a.handedness=='reversed':
        reflection=moon_plane_reflection(c,observed_moon)
        m=np.r_[reflection@m[:2],m[2]];s=np.r_[reflection@s[:2],s[2]]
        rotation=north_matrix[:,:2]@reflection
        north_matrix=np.column_stack([rotation,c-rotation@c])
    north_xy=np.stack([x,y],axis=-1)@north_matrix[:,:2].T+north_matrix[:,2]
    north_moon=north_matrix[:,:2]@observed_moon+north_matrix[:,2]
    moon=component(view,m,a.pmax);sun=component(view,s,a.pmax)
    with np.errstate(divide='ignore',invalid='ignore'):
        obs=np.hypot(Q,U)/I;phase=(Q+1j*U)/np.hypot(Q,U)
    disc=np.hypot(x-c[0],y-c[1])<.985*r
    mask=disc&np.isfinite(obs)&(I>0)&(obs>=0)&(obs<=.95)
    mask_hash=None
    if a.sky_mask:
        rawmask=np.load(a.sky_mask,allow_pickle=False)
        if rawmask.shape!=original_shape or rawmask.dtype!=np.bool_:p.error('Sky mask must be an original-shape boolean NPY')
        mask&=block_mean(rawmask.astype(float),a.bin_size)>.999
        mask_hash=digest(a.sky_mask)
    else:mask&=np.hypot(x-c[0],y-c[1])<site['eye_radius_px']
    for source in [m,s]:mask&=abs(np.sum(view*source,axis=-1))<np.cos(np.deg2rad(10))
    # Spatial blocks, not adjacent-pixel random splitting. Never optimize the split.
    held=((x//96+y//96).astype(int)%3)==0
    train=mask&~held;test=mask&held
    if min(train.sum(),test.sum())<30:p.error('Insufficient training or withheld sky samples')
    models={};arrays={}
    for kind in ['moon','sun','mixed']:
        result,intensity,z=fit_model(obs,moon,sun,train,kind)
        pred=abs(z)/intensity
        result.update(training=errors(obs,pred,train),withheld=errors(obs,pred,test))
        angle_mask=test&(obs>.05)&(abs(z)/intensity>.05)
        if min(I.shape)>4:
            predicted_phase=np.divide(z,abs(z),out=np.full(z.shape,np.nan+1j*np.nan),where=abs(z)>1e-12)
            result['withheld_aop_change_check_32px_equivalent']=residual_summary(phase,predicted_phase,angle_mask,max(1,32//a.bin_size))
        models[kind]=result;arrays[kind]=(intensity,z,pred)
    predicted_moon=pixels(m[None,:],c,r)[0]
    north_predicted_moon=north_matrix[:,:2]@predicted_moon+north_matrix[:,2]
    a.output.mkdir(parents=True)
    fig,axes=plt.subplots(2,4,figsize=(16,8),layout='constrained')
    panels=[('Measured camera DoLP',obs),('Lunar-only fitted dilution',arrays['moon'][2]),('Solar-only fitted dilution',arrays['sun'][2]),('Combined fitted dilution',arrays['mixed'][2])]
    lower=[('Training (blue) / withheld (red)',np.where(mask,np.where(held,1,-1),np.nan))]+[(k+' minus measured',np.where(mask,arrays[k][2]-obs,np.nan)) for k in ['moon','sun','mixed']]
    for row,items in enumerate([panels,lower]):
        for j,(title,data) in enumerate(items):
            ax=axes[row,j]
            im=ax.pcolormesh(north_xy[...,0],north_xy[...,1],np.where(disc,data,np.nan),shading='nearest',cmap='viridis' if row==0 else 'RdBu_r',vmin=0 if row==0 else (-1 if j==0 else -.25),vmax=1 if row==0 or j==0 else .25)
            ax.set_title(title,fontsize=10);ax.axis('off');fig.colorbar(im,ax=ax,shrink=.65)
            ax.add_patch(Circle(c,site['eye_radius_px'],fill=False,color='white',lw=.7))
            ax.set_aspect('equal');ax.set_xlim(c[0]-r,c[0]+r);ax.set_ylim(c[1]+r,c[1]-r)
            ax.plot(north_moon[0],north_moon[1],'kx',ms=6,mew=1.2)
            ax.plot(north_predicted_moon[0],north_predicted_moon[1],'+',color='red',ms=6,mew=1.2)
            ax.text(.03,.95,'N ↑',transform=ax.transAxes,fontsize=10)
            ax.text(.88,.5,'E →',transform=ax.transAxes,fontsize=8)
    fig.suptitle(a.site+'; handedness hypothesis: '+a.handedness+'; shared Moon-anchored north-up frame; fixed Sun/Moon geometry\nBlack ×: locked observed Moon; red +: astronomical Moon under provisional zenith-axis projection.\nOne common coordinate transform; no rotation fitted to polarization. Outside mask is extrapolation.',fontsize=11)
    fig.savefig(a.output/'twilight_moon_fit.png',dpi=150);plt.close(fig)
    np.savez_compressed(a.output/'twilight_moon_fit_arrays.npz',measured_dolp=obs,train=train,withheld=test,
        moon_dolp=arrays['moon'][2],sun_dolp=arrays['sun'][2],mixed_dolp=arrays['mixed'][2],mixed_I=arrays['mixed'][0],mixed_Q=arrays['mixed'][1].real,mixed_U=arrays['mixed'][1].imag)
    report=dict(handedness_hypothesis=a.handedness,coordinate_transform_determinant=float(np.linalg.det(north_matrix[:,:2])),reflection_hypothesis_applied=a.handedness=='reversed',display_frame='shared Moon-anchored north up; east right under stated parity hypothesis',display_camera_to_north_xy=north_matrix.tolist(),display_rotation_clockwise_deg=north_angle,observed_moon_north_xy=north_moon.tolist(),predicted_moon_north_xy=north_predicted_moon.tolist(),moon_pixel_separation=float(np.linalg.norm(observed_moon-predicted_moon)),site=a.site,status='exploratory_fit_not_production_model',models=models,
        parameter_definitions=dict(moon=['unpolarized_background / lunar_zenith_I'],sun=['unpolarized_background / solar_zenith_I'],mixed=['solar_zenith_I / lunar_zenith_I','unpolarized_background / lunar_zenith_I']),
        fixed_pmax=a.pmax,solar=solar,astronomical_anchor=site['astronomical_anchor'],
        source_lock_sha256=digest(a.reference_lock),stokes_sha256=digest(npz),script_sha256=digest(Path(__file__)),sky_mask_sha256=mask_hash,
        sky_mask='explicit user mask' if a.sky_mask else 'central 60 degree full field; no automatic sky segmentation',
        limitations=['Ideal Rayleigh-shaped components, not calibrated multiple-scattering twilight radiative transfer.',
        'Not the KS91/spectral production model; pmax is fixed and must be examined in sensitivity runs.',
        'Uniform unpolarized background is phenomenological dilution, not a fitted HPS/LED source or calibrated ALAN=5.',
        'Unknown instrument polarization and analyzer intensity inconsistency can bias the inferred mixture.',
        'Registration, camera handedness, optical-axis assumption and Moon identity remain provisional.',
        'Spatial withheld regions test interpolation in this capture; they are not independent nights.',
        'AoP changes are an unfitted diagnostic, not an absolute AoP calibration. No global threshold is estimated.',
        'Fit parameters may be non-identifiable; condition numbers and boundary flags are diagnostics, not confidence intervals.'])
    (a.output/'twilight_moon_fit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v['withheld'] for k,v in models.items()},indent=2));print('Saved',a.output)

if __name__=='__main__':main()
