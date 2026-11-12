"""Equation-traceable Wang et al. 2016 sky kernel, DOI 10.1088/2040-8978/18/6/065601.

Implements printed eqs 1,3-8,12-17 as numerical functions. Wavelength/basis
inputs are explicit: no fabricated supplementary coefficients. Complex w is
in the paper's stereographic coordinate convention; absolute AoP/camera basis
must be verified against measured skies before fitting. Solar-to-lunar use is
an explicit extension, not validation claimed by this solar paper.
"""
import numpy as np


def neutral_points(source_zenith,source_azimuth,separation_deg):
    """Printed eqs 4/5; delta=4 atan(L). Angles in radians except separation."""
    if not 0<separation_deg<90:raise ValueError('Explicit neutral-point separation must be between 0 and 90deg')
    L=np.tan(np.deg2rad(separation_deg)/4)
    r=np.tan(source_zenith/2);xs=r*np.cos(source_azimuth);ys=r*np.sin(source_azimuth)
    c=np.cos(source_azimuth);s=np.sin(source_azimuth)
    plus=(xs+L*c)/(1-L*c*xs)+1j*(ys+L*s)/(1-L*s*ys)
    minus=(xs-L*c)/(1+L*c*xs)+1j*(ys-L*s)/(1+L*s*ys)
    return complex(plus),complex(minus)


def complex_pattern(theta,phi,source_zenith,source_azimuth,separation_deg):
    """Eq3 on dimensionless stereographic z=tan(theta/2) exp(i phi)."""
    z=np.tan(theta/2)*np.exp(1j*phi)
    a,b=neutral_points(source_zenith,source_azimuth,separation_deg)
    if abs(a)<1e-12 or abs(b)<1e-12:
        raise ValueError('Printed eq3 singular at zero neutral coordinate; needs analytic-limit treatment')
    return ((z-a)*(z-b)*(z+1/a.conjugate())*(z+1/b.conjugate()) /
        ((1+abs(z)**2)**2 * abs(a+1/a.conjugate())*abs(b+1/b.conjugate())))


def turbidity_multiplier(turbidity,k1=4.,k2=.12):
    """Eq7: exp(-T/k1 + k2), not exp((-T+k2)/k1)."""
    if k1<=0 or turbidity<0:raise ValueError('Require k1>0 and T>=0')
    return float(np.exp(-turbidity/k1+k2))


def horizon_factor(theta,n=10.):
    if n<=0:raise ValueError('N must be positive')
    return np.maximum(np.cos(theta),0)**(1/n)


def perez(theta,gamma,coefficients):
    """Eq12; caller supplies verified A,B,C,D,E for luminance/x/y separately."""
    a,b,c,d,e=coefficients
    return (1+a*np.exp(b/np.maximum(np.cos(theta),1e-8)))*(1+c*np.exp(d*gamma)+e*np.cos(gamma)**2)


def cie_spectrum(x,y,s0,s1,s2):
    """Eqs13-15; spectral basis values must come from cited numerical tables."""
    den=.0241+.2562*x-.7341*y
    m1=(-1.3515-1.7703*x+5.9114*y)/den
    m2=(.0300-31.4424*x+30.0717*y)/den
    return s0+m1*s1+m2*s2


def spectral_factor(sd,s90,ssource):
    """Eq16; all spectra need a common normalization and wavelength grid."""
    if np.any(np.asarray(sd)<=0) or np.any(np.asarray(ssource)<=0) or np.any(np.asarray(ssource)==s90):
        raise ValueError('Eq16 has nonpositive or degenerate spectral inputs')
    return (1/np.asarray(sd)-1/np.asarray(ssource))*np.asarray(s90)*ssource/(ssource-s90)


def natural_stokes(theta,phi,source_zenith,source_azimuth,separation_deg,
                   intensity,turbidity=3.,mode='eq6_horizon',spectral_s=None,
                   depolarization_multiplier=None):
    w=complex_pattern(theta,phi,source_zenith,source_azimuth,separation_deg)
    if mode=='eq17':
        if spectral_s is None:raise ValueError('Full eq17 requires explicitly supplied spectral S from eq16')
        weight=theta*horizon_factor(theta)+(np.pi/2-theta)*np.asarray(spectral_s)
    elif mode=='eq6_horizon':
        # Explicit reduced pilot: eq6 multiplied by eq8, NOT full eq17.
        weight=horizon_factor(theta)
    else:raise ValueError('Choose eq17 or explicit reduced eq6_horizon pilot')
    m=turbidity_multiplier(turbidity) if depolarization_multiplier is None else float(depolarization_multiplier)
    if m<0:raise ValueError('Depolarization multiplier cannot be negative')
    i=np.broadcast_to(np.asarray(intensity,dtype=float),np.shape(theta)).copy()
    p=abs(w)*weight*m
    if np.any(~np.isfinite(p)) or np.any(p<0) or np.any(p>1+1e-10) or np.any(i<=0):
        raise ValueError('Nonphysical natural Stokes parameters; inspect paper inputs/conventions, do not clip silently')
    angle=.5*np.angle(w)
    return i,i*p*np.cos(2*angle),i*p*np.sin(2*angle)


def add_unpolarized(i,q,u,artificial_intensity):
    addition=np.asarray(artificial_intensity,dtype=float)
    if np.any(~np.isfinite(addition)) or np.any(addition<0):raise ValueError('ALAN must be finite and nonnegative')
    return np.asarray(i)+addition,np.asarray(q).copy(),np.asarray(u).copy()


def equivalent_turbidity_increment(relative_addition,k1=4.):
    """Only when I_ALAN/I_natural is spatially constant; NOT an atmospheric identity."""
    if relative_addition<0:raise ValueError('Relative addition must be nonnegative')
    return float(k1*np.log1p(relative_addition))


def lunar_track(utc_times,latitude,longitude,height_m=0.):
    """Astropy built-in ephemeris, explicit UTC, pressure=0 (no refraction).

    Does not retrieve new IERS data. Timestamps must be independently verified.
    Moon intensity/phase/spectrum not derived by this geometric function.
    """
    from astropy import units as units
    from astropy.coordinates import EarthLocation,AltAz,get_body,solar_system_ephemeris
    from astropy.time import Time
    from astropy.utils import iers
    iers.conf.auto_download=False
    location=EarthLocation(lat=latitude*units.deg,lon=longitude*units.deg,height=height_m*units.m)
    times=Time(utc_times,scale='utc')
    with solar_system_ephemeris.set('builtin'):
        coords=get_body('moon',times,location).transform_to(AltAz(obstime=times,location=location,pressure=0*units.hPa))
    return coords.alt.deg,coords.az.deg
