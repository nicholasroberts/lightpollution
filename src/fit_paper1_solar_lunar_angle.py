"""Fit a fixed-geometry Berry/Wang Sun/Moon mixture to spatial AoP changes.

No ALAN, DoLP ceiling, global analyzer-zero, sky rotation or parity fit.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model_paper1_lunar_baseline import lunar_component
from audit_spider_lunar_aop import block_mean
from spider_fisheye_registration import rays
from rotate_spider_measured_north import digest


def optimize_fraction(measured_change, lunar_lo, lunar_hi, solar_lo, solar_hi, selected):
    def residual(f):
        lo=(1-f)*lunar_lo+f*solar_lo
        hi=(1-f)*lunar_hi+f*solar_hi
        delta=.5*np.angle(measured_change*np.conj(hi*np.conj(lo)))
        ambiguous=(abs(lo)<1e-10)|(abs(hi)<1e-10)
        return np.where(ambiguous,np.pi/2,delta)
    def cost(f):
        r=residual(f)[selected]
        return float(np.mean(r*r))
    grid=np.linspace(0,1,401)
    costs=np.array([cost(f) for f in grid])
    minima=[j for j in range(1,len(grid)-1) if costs[j]<=costs[j-1] and costs[j]<=costs[j+1]]
    candidates=[(costs[0],0.),(costs[-1],1.)]
    for j in minima:
        result=minimize_scalar(cost,bounds=(grid[j-1],grid[j+1]),method='bounded',
                               options={'xatol':1e-10})
        candidates.append((result.fun,float(result.x)))
    best=min(candidates)
    return best[1],residual,grid,costs


def error_stats(values):
    v=np.rad2deg(values)
    return dict(pair_count=int(v.size),median_absolute_deg=float(np.median(abs(v))),
        rms_deg=float(np.sqrt(np.mean(v*v))),p95_absolute_deg=float(np.percentile(abs(v),95)),
        mean_signed_deg=float(v.mean()))


def sun_at_capture(anchor,site,reference):
    from astropy.coordinates import EarthLocation,AltAz,get_sun
    from astropy.time import Time
    from astropy.utils import iers
    import astropy.units as units
    iers.conf.auto_download=False
    rows=[r for r in csv.DictReader(reference.open()) if Path(r['source_filename']).stem==site]
    if len(rows)!=1:raise ValueError('Need exactly one verified site row')
    row=rows[0]
    if anchor.get('site_id') and anchor['site_id']!=row['site_id']:
        raise ValueError('Site ID and astronomical anchor differ')
    utc=anchor['utc']
    moment=Time(utc.replace('Z',''),scale='utc')
    lat,lon=float(row['latitude']),float(row['longitude'])
    location=EarthLocation(lat=lat*units.deg,lon=lon*units.deg,height=0*units.m)
    sun=get_sun(moment).transform_to(AltAz(obstime=moment,location=location,pressure=0*units.hPa))
    return dict(sun_elevation_deg=float(sun.alt.deg),sun_compass_azimuth_deg=float(sun.az.deg),
        utc=utc,latitude_deg=lat,longitude_deg=lon,height_m=0,pressure_hpa=0,
        source='astropy.get_sun; offline IERS; geometric AltAz')


def source_vector(elevation,azimuth):
    el,az=np.deg2rad([elevation,azimuth])
    return np.array([-np.cos(el)*np.sin(az),np.cos(el)*np.cos(az),np.sin(el)])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-lock',type=Path,required=True)
    p.add_argument('--north',type=Path,required=True)
    p.add_argument('--site',default='museum')
    p.add_argument('--reference',type=Path,default=Path(__file__).resolve().parents[1]/'data/reference/prc_spider_18_sites_verified.csv')
    p.add_argument('--neutral-separation-deg',type=float,required=True)
    p.add_argument('--bin-size',type=int,default=8)
    p.add_argument('--lags-px',type=int,nargs='+',default=[32,64])
    p.add_argument('--minimum-dolp',type=float,default=.05)
    p.add_argument('--block-px',type=int,default=96)
    p.add_argument('--sky-mask',type=Path,help='Optional original north-array shape boolean clear-sky mask')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use fresh output directory')
    if a.bin_size<1 or a.block_px<a.bin_size or not 0<a.minimum_dolp<1 or not 0<a.neutral_separation_deg<90:
        p.error('Invalid scale, threshold or neutral separation')
    if any(l<a.bin_size or l%a.bin_size for l in a.lags_px):
        p.error('Lags must be positive integer multiples of bin size')
    nrpath=a.north/'north_rotation.json';npz=a.north/'north_stokes.npz'
    nr=json.loads(nrpath.read_text())
    if nr['site']!=a.site or nr['source_lock_sha256']!=digest(a.reference_lock):
        p.error('North product and source lock differ')
    lock=json.loads(a.reference_lock.read_text())['sites'][a.site]
    with np.load(npz,allow_pickle=False) as z:
        original_shape=z['I'].shape
        intensity,q,u=[block_mean(z[k],a.bin_size) for k in ['I','Q','U']]
        center=np.array(z['center_xy'])
    y,x=np.indices(intensity.shape,dtype=float)
    x=x*a.bin_size+(a.bin_size-1)/2;y=y*a.bin_size+(a.bin_size-1)/2
    view=rays(np.stack([x,2*center[1]-y],-1),center,lock['rim_radius_px'])
    anchor=lock['astronomical_anchor'];solar=sun_at_capture(anchor,a.site,a.reference)
    moon=source_vector(anchor['moon_elevation_deg'],anchor['moon_compass_azimuth_deg'])
    sun=source_vector(solar['sun_elevation_deg'],solar['sun_compass_azimuth_deg'])
    mi,mz=lunar_component(view,moon,a.neutral_separation_deg)
    si,sz=lunar_component(view,sun,a.neutral_separation_deg)
    observed=q+1j*u
    with np.errstate(divide='ignore',invalid='ignore'):
        dolp=abs(observed)/intensity
    mask=(np.hypot(x-center[0],y-center[1])<lock['eye_radius_px'])&(intensity>0)&np.isfinite(observed)&np.isfinite(dolp)&(dolp>=a.minimum_dolp)&(dolp<=.95)
    # Fixed support for endpoints and every candidate fraction; no model-specific threshold masks.
    mask&=np.isfinite(mz)&np.isfinite(sz)&(abs(mz)>1e-8)&(abs(sz)>1e-8)
    if a.sky_mask:
        raw=np.load(a.sky_mask,allow_pickle=False)
        if raw.shape!=original_shape or raw.dtype!=np.bool_:p.error('Sky mask shape/dtype differs')
        mask&=block_mean(raw.astype(float),a.bin_size)>.999
    for s in [moon,sun]:
        mask&=abs(np.sum(view*s,axis=-1))<np.cos(np.deg2rad(10))
    bx=(x//a.block_px).astype(int);by=(y//a.block_px).astype(int)
    withheld=(bx+by)%3==0
    block=bx+by*(int(bx.max())+1)
    pairs=[]
    for spacing in sorted(set(a.lags_px)):
        lag=spacing//a.bin_size
        if lag>=min(mask.shape):p.error('Lag too large')
        for axis in [0,1]:
            lo=[slice(None)]*2;hi=lo.copy();lo[axis]=slice(None,-lag);hi[axis]=slice(lag,None)
            lo,hi=tuple(lo),tuple(hi)
            good=mask[lo]&mask[hi]
            # Both endpoints must be in the same training/withheld partition.
            good&=withheld[lo]==withheld[hi]
            pairs.append((observed[hi][good]*np.conj(observed[lo][good]),
                mz[lo][good],mz[hi][good],sz[lo][good],sz[hi][good],
                withheld[lo][good],block[lo][good],block[hi][good],
                np.full(int(good.sum()),spacing)))
    columns=[np.concatenate([p[j] for p in pairs]) for j in range(9)]
    measured,ml,mh,sl,sh,test,bl,bh,spacings=columns
    train=~test
    if min(train.sum(),test.sum())<30:p.error('Insufficient training or withheld pairs')
    fraction,residual,grid,costs=optimize_fraction(measured,ml,mh,sl,sh,train)
    zi=(1-fraction)*mz+fraction*sz
    ii=(1-fraction)*mi+fraction*si
    metrics={}
    for name,f in [('lunar_only',0.),('solar_only',1.),('combined',fraction)]:
        r=residual(f)
        metrics[name]={part:{str(s):error_stats(r[selection&(spacings==s)]) for s in sorted(set(spacings.tolist()))
                            if (selection&(spacings==s)).any()}
                        for part,selection in [('training',train),('withheld',test)]}
        metrics[name]['pooled_training']=error_stats(r[train])
        metrics[name]['pooled_withheld']=error_stats(r[test])
    blocks=np.unique(np.r_[bl[train],bh[train]])
    stability=[]
    for bid in blocks:
        keep=train&(bl!=bid)&(bh!=bid)
        if keep.sum()<30:continue
        f,_,_,_=optimize_fraction(measured,ml,mh,sl,sh,keep)
        stability.append(dict(omitted_block=int(bid),solar_fraction=f))
    a.output.mkdir(parents=True)
    with (a.output/'ratio_profile.csv').open('w',newline='') as stream:
        w=csv.writer(stream);w.writerow(['solar_fraction','solar_to_lunar_ratio','training_AoP_change_MSE_rad2'])
        for f,cost in zip(grid,costs):w.writerow([f,f/(1-f) if f<1 else 'solar_only',cost])
    np.savez_compressed(a.output/'solar_lunar_fit_arrays.npz',measured_I=intensity,measured_Q=q,measured_U=u,
        lunar_I=mi,lunar_Q=mz.real,lunar_U=mz.imag,solar_I=si,solar_Q=sz.real,solar_U=sz.imag,
        combined_I=ii,combined_Q=zi.real,combined_U=zi.imag,
        fit_mask=mask,withheld_pixel_partition=withheld,x_px=x,y_px=y)
    fig,axes=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
    axes[0].plot(grid,np.rad2deg(np.sqrt(costs)),color='#0072B2')
    axes[0].axvline(fraction,color='black',ls='--')
    axes[0].set(xlabel='Solar fraction of uniform natural intensity',ylabel='Training RMS AoP-change error (degrees)',title='Mixture profile; no analyzer-zero fit')
    names=list(metrics);axes[1].bar(names,[metrics[n]['pooled_withheld']['rms_deg'] for n in names],color=['#777777','#E69F00','#0072B2'])
    axes[1].set(ylabel='Withheld RMS AoP-change error (degrees)',title='Identical held-out pairs')
    fig.suptitle(a.site+' — solar/lunar fit diagnostics; not the assembled publication panel')
    fig.savefig(a.output/'fit_diagnostics.png',dpi=180);plt.close(fig)
    result=dict(site=a.site,status='angle_structure_fit_no_ALAN',solar_fraction=fraction,
        lunar_fraction=1-fraction,solar_to_lunar_uniform_intensity_ratio=fraction/(1-fraction) if fraction<1 else None,
        normalization='Each component has uniform I=1, as in the current Berry/Wang baseline. Fraction is camera-green proxy intensity weight, not top-of-atmosphere or electrical source power.',
        objective='Mean squared circular AoP-change residual; horizontal and vertical pairs, fixed mask, training blocks only. Constant analyzer zero cancels.',
        source_lock_sha256=digest(a.reference_lock),source_north_arrays_sha256=digest(npz),
        source_north_report_sha256=digest(nrpath),script_sha256=digest(__file__),
        kernel_sha256=digest(Path(__file__).with_name('wang2016_sky.py')),
        component_module_sha256=digest(Path(__file__).with_name('model_paper1_lunar_baseline.py')),
        site_reference_sha256=digest(a.reference),
        sky_mask_sha256=digest(a.sky_mask) if a.sky_mask else None,
        astronomical_anchor=anchor,solar=solar,neutral_separation_deg=a.neutral_separation_deg,
        bin_size=a.bin_size,lags_px=sorted(set(a.lags_px)),block_px=a.block_px,minimum_measured_dolp=a.minimum_dolp,
        central_mask_pixels=int(mask.sum()),training_pair_count=int(train.sum()),withheld_pair_count=int(test.sum()),
        training_blocks=int(len(blocks)),metrics=metrics,leave_one_training_block_out=stability,
        fraction_on_boundary=bool(fraction<1e-6 or fraction>1-1e-6),
        fitted_background=0,fitted_pattern_rotation=False,fitted_analyzer_zero=False,
        limitations=['Berry/Wang geometric solar component extended below the horizon; not calibrated twilight radiative transfer.',
        'No separate solar/lunar spectra: fraction is a green-channel effective weight under achromatic component assumptions.',
        'Correlated pair errors and overlapping lags are not independent replicates; no naive pixel p-values or confidence intervals.',
        'Leave-one-block-out results measure conditional spatial stability, not a formal confidence interval or between-night uncertainty.',
        'Mask is central physical sky candidate unless explicit clear-sky mask supplied; no automatic foreground segmentation.',
        'Model selection on this one capture is exploratory; withheld blocks test spatial interpolation, not independent nights.',
        'No ALAN or amplitude fit; combined DoLP remains a prediction and is not optimized.',
        'Inherited camera parity, unknown analyzer zero and idealized optical-axis/projection remain uncertainties.'])
    (a.output/'solar_lunar_fit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(solar=solar,solar_fraction=fraction,solar_to_lunar_ratio=result['solar_to_lunar_uniform_intensity_ratio'],
        withheld_metrics={n:metrics[n]['pooled_withheld'] for n in metrics}),indent=2))
    print('Saved',a.output)


if __name__=='__main__':main()
