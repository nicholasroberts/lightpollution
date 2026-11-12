"""Nightly fixed-body-field PRC scenarios; not fitted ALAN or absolute lunar brightness.

Uses Berry-normalized geometry + Wang horizon factor. Spectral input columns:
wavelength_nm,natural,alan,sensitivity; optional dop_scale,abs_parallel,abs_perpendicular.
Natural/ALAN spectra normalized to unit unweighted band integral separately.
ALAN ratios are radiometric band-integral ratios, not VIIRS radiance or lux.
"""
from pathlib import Path
from datetime import datetime,timedelta,timezone
from zoneinfo import ZoneInfo
import argparse,csv,json
import numpy as np
from scipy.integrate import trapezoid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wang2016_sky import berry_normalized_pattern,horizon_factor
from reconstruct_spider_eye_contrast import digest
ROOT=Path(__file__).resolve().parents[1]


def spectral_coefficients(path=None):
    if path is None:
        r,t=10**(-.4),10**(-.04)
        return (1-(r+t)/2,(t-r)/2,1-(r+t)/2),{'mode':'ACHROMATIC_SCENARIO_NOT_SPECIES_SPECTRAL_MODEL'}
    data=np.genfromtxt(path,delimiter=',',names=True)
    if data.ndim!=1 or len(data)<2:raise ValueError('Need at least two spectral samples')
    names=data.dtype.names
    for name in ['wavelength_nm','natural','alan','sensitivity']:
        if name not in names:raise ValueError('Missing spectral column '+name)
    wave=data['wavelength_nm']
    if np.any(~np.isfinite(wave)) or np.any(np.diff(wave)<=0):raise ValueError('Wavelengths must increase')
    get=lambda name,default: data[name] if name in names else np.full(wave.shape,default)
    n,a,s=[data[k] for k in ['natural','alan','sensitivity']]
    d=get('dop_scale',1.);ap=get('abs_parallel',.4);at=get('abs_perpendicular',.04)
    if any(np.any(~np.isfinite(v)) or np.any(v<0) for v in [n,a,s,d,ap,at]) or np.any(d>1) or np.any(at>ap):raise ValueError('Invalid nonnegative spectra, DoP scale or absorbances')
    integ=lambda x:float(trapezoid(x,wave))
    if min(integ(n),integ(a),integ(s))<=0:raise ValueError('Spectra and sensitivity must have positive integrals')
    n=n/integ(n);a=a/integ(a)
    b=1-(10**(-ap)+10**(-at))/2;g=(10**(-at)-10**(-ap))/2
    coefficients=(integ(s*n*b),integ(s*n*g*d),integ(s*a*b))
    if min(coefficients[0],coefficients[2])<=0:raise ValueError('Nonpositive receptor response')
    return coefficients,{'mode':'SUPPLIED_SPECTRAL_RESPONSE','path':str(path),'sha256':digest(path),'wavelength_range_nm':[float(wave[0]),float(wave[-1])],'samples':len(wave)}


def body_field(heading_deg,tilt_deg,samples=1024):
    """Uniform-solid-angle cap; tilt from zenith towards body heading (north=0)."""
    h,t=np.deg2rad([heading_deg,tilt_deg]);forward=np.array([np.sin(h),np.cos(h),0.]);up=np.array([0.,0.,1.])
    center=np.cos(t)*up+np.sin(t)*forward
    axis=np.cos(t)*forward-np.sin(t)*up;side=np.cross(center,axis)
    k=np.arange(samples);cosrho=1-(k+.5)/samples*(1-np.cos(np.pi/6));sinrho=np.sqrt(1-cosrho**2)
    angle=k*np.pi*(3-np.sqrt(5))
    rays=cosrho[:,None]*center+sinrho[:,None]*(np.cos(angle)[:,None]*axis+np.sin(angle)[:,None]*side)
    # Parallel transport body analyzer from cap centre to each viewing ray.
    transported=axis-(rays@axis/(1+rays@center))[:,None]*(center+rays)
    return rays,transported


def integrated_pattern(rays,analyzer,moon_elevation,moon_azimuth,separation):
    if np.any(rays[:,2]<-1e-12):raise ValueError('FOV extends below horizon; choose tilt <=60deg')
    theta=np.arccos(np.clip(rays[:,2],-1,1));phi=np.arctan2(rays[:,1],rays[:,0])
    # Model phi is mathematical azimuth from east towards north.
    w=berry_normalized_pattern(theta,phi,np.deg2rad(90-moon_elevation),np.deg2rad(90-moon_azimuth),separation)
    w=w*horizon_factor(theta)
    et=np.column_stack([np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),-np.sin(theta)])
    ep=np.column_stack([-np.sin(phi),np.cos(phi),np.zeros_like(phi)])
    ex=np.cos(phi)[:,None]*et-np.sin(phi)[:,None]*ep
    ey=np.sin(phi)[:,None]*et+np.cos(phi)[:,None]*ep
    z=np.sum(analyzer*ex,axis=1)+1j*np.sum(analyzer*ey,axis=1)
    return float(np.mean(np.real(w*np.conj(z)**2)))


def response_contrast(pattern,atmospheric_factor,alan_ratio,coefficients):
    bn,g,ba=coefficients
    mean=bn+alan_ratio*ba;mod=g*atmospheric_factor*pattern
    if mean<=abs(mod):raise ValueError('Nonphysical receptor response')
    return float(np.log((mean-mod)/(mean+mod)))


def night_track(anchor,latitude,longitude,step):
    from astropy.coordinates import EarthLocation,AltAz,get_body,get_sun,solar_system_ephemeris
    from astropy.time import Time
    from astropy import units as u
    from astropy.utils import iers
    iers.conf.auto_download=False
    zone=ZoneInfo(anchor['timezone']);capture=datetime.fromisoformat(anchor['utc']).replace(tzinfo=timezone.utc)
    local=capture.astimezone(zone);noon=local.replace(hour=12,minute=0,second=0,microsecond=0)
    if local.hour<12:noon-=timedelta(days=1)
    start=noon.astimezone(timezone.utc);end=(noon+timedelta(days=1)).astimezone(timezone.utc)
    stamps=[start+timedelta(minutes=k) for k in range(0,int((end-start).total_seconds()/60),step)]
    times=Time(stamps);location=EarthLocation(lat=latitude*u.deg,lon=longitude*u.deg)
    frame=AltAz(obstime=times,location=location,pressure=0*u.hPa)
    with solar_system_ephemeris.set('builtin'):moon=get_body('moon',times,location).transform_to(frame)
    sun=get_sun(times).transform_to(frame)
    return stamps,moon.alt.deg,moon.az.deg,sun.alt.deg


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--moon-anchors',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--prc-colour-limit',type=float,default=2.,help='Shared symmetric PRC colour limit for every site and atmosphere (default +/-2)')
    p.add_argument('--spectra',type=Path);p.add_argument('--references',type=Path,default=ROOT/'data/reference/prc_spider_18_sites_verified.csv')
    p.add_argument('--atmospheric-factors',type=float,nargs='+',default=[1.,.8,.6])
    p.add_argument('--alan-ratios',type=float,nargs='+',default=[0.,.5,1.,3.])
    p.add_argument('--tilts',type=float,nargs='+',default=[0.,30.,60.]);p.add_argument('--step-minutes',type=int,default=20)
    p.add_argument('--heading-step',type=int,default=10);p.add_argument('--samples',type=int,default=1024)
    p.add_argument('--neutral-separation-deg',type=float,default=24.)
    args=p.parse_args()
    if args.output.exists():p.error('Use fresh output directory')
    if not np.isfinite(args.prc_colour_limit) or args.prc_colour_limit<=0:p.error('Colour limit must be positive and finite')
    if not 0<args.step_minutes<=120 or not 0<args.heading_step<=90 or args.samples<64:p.error('Invalid sampling')
    if any(not np.isfinite(x) or not 0<=x<=1 for x in args.atmospheric_factors):p.error('Atmospheric factors must be in [0,1]')
    if any(not np.isfinite(x) or x<0 for x in args.alan_ratios):p.error('ALAN ratios must be nonnegative')
    if any(not np.isfinite(x) or not 0<=x<=60 for x in args.tilts):p.error('Tilts must be 0..60deg')
    coefficients,spectral_info=spectral_coefficients(args.spectra)
    anchors=json.loads(args.moon_anchors.read_text())['sites'];refs={r['site_id']:r for r in csv.DictReader(args.references.open())}
    headings=np.arange(0,360,args.heading_step);rows=[];summaries=[];tracks=[]
    args.output.mkdir(parents=True)
    for site in ['Cottage','oxfordbussinesspark']:
        anchor=anchors[site];ref=refs[anchor['site_id']]
        stamps,elevations,azimuths,sun=night_track(anchor,float(ref['latitude']),float(ref['longitude']),args.step_minutes)
        zone=ZoneInfo(anchor['timezone']);hours=np.array([(t-stamps[0]).total_seconds()/3600+12 for t in stamps])
        print(f'{site}: computing {len(stamps)} times, {len(headings)} headings',flush=True)
        matrices={}
        for tilt in args.tilts:
            base=np.full((len(headings),len(stamps)),np.nan)
            for ih,heading in enumerate(headings):
                rays,analyzer=body_field(heading,tilt,args.samples)
                for it,(stamp,el,az,solar) in enumerate(zip(stamps,elevations,azimuths,sun)):
                    if solar>=-6 or not .1<el<89.9:continue
                    base[ih,it]=integrated_pattern(rays,analyzer,el,az,args.neutral_separation_deg)
            for atmosphere in args.atmospheric_factors:
                for alan in args.alan_ratios:
                    matrix=np.full_like(base,np.nan)
                    for ih,it in zip(*np.where(np.isfinite(base))):
                        value=response_contrast(base[ih,it],atmosphere,alan,coefficients);matrix[ih,it]=value
                        rows.append(dict(site=site,utc=stamps[it].isoformat(),local_time=stamps[it].astimezone(zone).isoformat(),moon_elevation_deg=float(elevations[it]),moon_azimuth_deg=float(azimuths[it]),body_heading_deg=int(headings[ih]),body_tilt_deg=tilt,atmospheric_factor=atmosphere,alan_band_ratio=alan,signed_prc=value,absolute_prc=abs(value)))
                    matrices[tilt,atmosphere,alan]=matrix
                    for it in range(len(stamps)):
                        vals=abs(matrix[:,it]);vals=vals[np.isfinite(vals)]
                        if len(vals):summaries.append(dict(site=site,utc=stamps[it].isoformat(),body_tilt_deg=tilt,atmospheric_factor=atmosphere,alan_band_ratio=alan,min_abs_prc=float(vals.min()),median_abs_prc=float(np.median(vals)),max_abs_prc=float(vals.max())))
        # One figure per atmospheric factor, all posture/ALAN combinations.
        for atmosphere in args.atmospheric_factors:
            fig,axes=plt.subplots(len(args.tilts),len(args.alan_ratios),figsize=(4*len(args.alan_ratios),3*len(args.tilts)),squeeze=False)
            limit=args.prc_colour_limit
            if any(np.any(abs(m[np.isfinite(m)])>limit) for key,m in matrices.items() if key[1]==atmosphere):
                raise ValueError('PRC exceeds shared display limit; increase --prc-colour-limit for all figures')
            for ri,tilt in enumerate(args.tilts):
                for ci,alan in enumerate(args.alan_ratios):
                    ax=axes[ri,ci];cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('lightgray')
                    art=ax.pcolormesh(hours,headings,matrices[tilt,atmosphere,alan],shading='nearest',cmap=cm,vmin=-limit,vmax=limit)
                    ax.set(title=f'Tilt {tilt:g}deg; ALAN/natural {alan:g}',xlabel='Local hour (24 = midnight)',ylabel='Body heading from north (deg)')
            night_hours=hours[sun<-6]
            if len(night_hours):
                for ax in axes.ravel():ax.set_xlim(night_hours.min()-.25,night_hours.max()+.25)
            fig.subplots_adjust(top=.88,hspace=.5,wspace=.45,right=.86)
            color_axis=fig.add_axes([.89,.2,.015,.6])
            fig.colorbar(art,cax=color_axis,label='Signed eye log contrast')
            fig.suptitle(f'{site}: atmospheric factor {atmosphere:g}; '+spectral_info['mode']+'\nGray = Moon below horizon / Sun above -6deg; fixed body field 60deg full width',fontsize=10)
            
            fig.savefig(args.output/f'{site}_atmosphere_{atmosphere:g}_heading.png',dpi=140);plt.close(fig)
        fig,axes=plt.subplots(len(args.tilts),1,figsize=(11,3*len(args.tilts)),squeeze=False)
        for ax,tilt in zip(axes[:,0],args.tilts):
            for alan in args.alan_ratios:
                m=abs(matrices[tilt,args.atmospheric_factors[0],alan]);valid=np.isfinite(m).any(axis=0)
                median=np.full(len(stamps),np.nan);lo=median.copy();hi=median.copy()
                if valid.any():median[valid]=np.median(m[:,valid],axis=0);lo[valid]=np.min(m[:,valid],axis=0);hi[valid]=np.max(m[:,valid],axis=0)
                line,=ax.plot(hours,median,label=f'ALAN {alan:g}');ax.fill_between(hours,lo,hi,color=line.get_color(),alpha=.12)
            ax.set(title=f'Tilt {tilt:g}deg; atmospheric factor {args.atmospheric_factors[0]:g}',xlabel='Local hour (24 = midnight)',ylabel='Absolute PRC');ax.legend()
        fig.suptitle(f'{site}: median and min–max over sampled headings (scenario, not confidence interval)');fig.tight_layout(rect=[0,0,1,.95]);fig.savefig(args.output/f'{site}_prc_ranges.png',dpi=140);plt.close(fig)
        tracks.append(dict(site=site,anchor=anchor,reference=ref,valid_moonlit_night_samples=int(np.sum((sun<-6)&(elevations>.1)&(elevations<89.9)))))
    for name,records in [('heading_responses.csv',rows),('prc_ranges.csv',summaries)]:
        if records:
            with (args.output/name).open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=records[0]);writer.writeheader();writer.writerows(records)
    payload=dict(status='FORWARD_SCENARIOS_NOT_FITTED_SITE_ALAN',parameters={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},spectra=spectral_info,sites=tracks,script_sha256=digest(__file__),limitations=['Natural radiance constant over directions and night; lunar phase/extinction/brightness not simulated','ALAN spatially uniform absolute radiance and constant over night; its ratio uses constant natural baseline','Atmospheric factor scales Q/U only; sensitivity analysis, not measured turbidity','No terrain/cloud masks or twilight polarization; modeled Moon above horizon and Sun below -6deg only','Body tilts and eye alignment are hypotheses; analyzer axes parallel transported over field','Berry unity-peak intensity interpreted as DoLP envelope; no restoration of Wang spectral eq17','Spectra separable from spatial pattern, fixed across night; wavelength effects only as supplied','No behavioral threshold or probabilities assigned to headings; ranges describe sampled scenarios','Capture dates and site matching retain input flags; Cottage date discrepancy unresolved'])
    (args.output/'provenance.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    print(f'Wrote {len(rows)} heading responses to {args.output}',flush=True)

if __name__=='__main__':main()
