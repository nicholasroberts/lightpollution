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


# Wang supplement S1-S3 (D_Y intercept retained exactly as Wang: -2.5571).
PEREZ_MATRICES={
 'Y':np.array([[.1787,-1.4630],[-.3554,.4275],[-.0227,5.3251],[.1206,-2.5571],[-.0670,.3703]]),
 'x':np.array([[.0193,-.2592],[-.0665,.0008],[-.0004,.2125],[-.0641,-.8989],[-.0033,.0452]]),
 'y':np.array([[-.0167,-.2608],[-.0950,.0092],[-.0079,.2102],[-.0441,-1.6537],[-.0109,.0529]])}
# Corrected S4/S5: Preetham author manuscript Appendix A.2, as agreed with investigator.
ZENITH_X=np.array([[.00166,-.00375,.00209,0],[-.02903,.06377,-.03202,.00394],[.11693,-.21196,.06052,.25886]])
ZENITH_Y=np.array([[.00275,-.00610,.00317,0],[-.04214,.08970,-.04153,.00516],[.15346,-.26756,.06670,.26688]])
PRINTED_DUPLICATE=np.array([[.0028,-.0061,.0032,0],[-.0421,.0897,-.0415,.0052],[.1535,-.2676,.0667,.2669]])
# Preetham Table 2: relative daylight basis, 380..780nm at 10nm. No extrapolation.
DAYLIGHT_BASIS=np.array([
[63.4,38.5,3],[65.8,35,1.2],[94.8,43.4,-1.1],[104.8,46.3,-.5],
[105.9,43.9,-.7],[96.8,37.1,-1.2],[113.9,36.7,-2.6],[125.6,35.9,-2.9],
[125.5,32.6,-2.8],[121.3,27.9,-2.6],[121.3,24.3,-2.6],[113.5,20.1,-1.8],
[113.1,16.2,-1.5],[110.8,13.2,-1.3],[106.5,8.6,-1.2],[108.8,6.1,-1],
[105.3,4.2,-.5],[104.4,1.9,-.3],[100,0,0],[96,-1.6,.2],
[95.1,-3.5,.5],[89.1,-3.5,2.1],[90.5,-5.8,3.2],[90.3,-7.2,4.1],
[88.4,-8.6,4.7],[84,-9.5,5.1],[85.1,-10.9,6.7],[81.9,-10.7,7.3],
[82.6,-12,8.6],[84.9,-14,9.8],[81.3,-13.6,10.2],[71.9,-12,8.3],
[74.3,-13.3,9.6],[76.4,-12.9,8.5],[63.3,-10.6,7],[71.7,-11.6,7.6],
[77,-12.2,8],[65.2,-10.2,6.7],[47.7,-7.8,5.2],[68.6,-11.2,7.4],[65,-10.4,6.8]])


def zenith_chromaticity(source_zenith,turbidity,convention='corrected_preetham'):
    t=np.array([turbidity**2,turbidity,1.]);z=np.array([source_zenith**3,source_zenith**2,source_zenith,1.])
    if convention=='corrected_preetham':mx,my=ZENITH_X,ZENITH_Y
    elif convention=='supplement_as_printed':mx=my=PRINTED_DUPLICATE
    else:raise ValueError('Unknown chromaticity convention')
    return float(t@mx@z),float(t@my@z)


def scattering_angle(theta,phi,source_zenith,source_azimuth):
    return np.arccos(np.clip(np.cos(theta)*np.cos(source_zenith)+
        np.sin(theta)*np.sin(source_zenith)*np.cos(phi-source_azimuth),-1,1))


def sky_chromaticity(theta,gamma,source_zenith,turbidity,convention='corrected_preetham'):
    xz,yz=zenith_chromaticity(source_zenith,turbidity,convention)
    t=np.array([turbidity,1.]);xcoeff=PEREZ_MATRICES['x']@t;ycoeff=PEREZ_MATRICES['y']@t
    return (xz*perez(theta,gamma,xcoeff)/perez(0.,source_zenith,xcoeff),
            yz*perez(theta,gamma,ycoeff)/perez(0.,source_zenith,ycoeff))


def relative_daylight_spectrum(theta,gamma,source_zenith,turbidity,wavelength_nm,convention='corrected_preetham'):
    wavelength_nm=float(wavelength_nm)
    if not 380<=wavelength_nm<=780:raise ValueError('Verified basis covers 380..780nm; no extrapolation')
    basis=[np.interp(wavelength_nm,np.arange(380,781,10),DAYLIGHT_BASIS[:,k]) for k in range(3)]
    x,y=sky_chromaticity(theta,gamma,source_zenith,turbidity,convention)
    return cie_spectrum(x,y,*basis)


def wang_spectral_inputs(theta,phi,source_zenith,source_azimuth,turbidity,wavelength_nm,
                         convention='corrected_preetham'):
    """Eqs10-16 plus explicit 90deg reference direction in source-zenith meridian.

    Reference ray opposite source azimuth has theta90=pi/2-theta_source.
    This choice is explicit: 90deg scattering alone does not uniquely specify
    a sky direction when intensity/chromaticity also depend on zenith angle.
    Source spectrum here is circumsource SKYLIGHT chromaticity, not direct disc.
    """
    if not 0<source_zenith<np.pi/2:raise ValueError('Spectral reference requires source strictly above horizon and below zenith')
    gamma=scattering_angle(theta,phi,source_zenith,source_azimuth)
    sd=relative_daylight_spectrum(theta,gamma,source_zenith,turbidity,wavelength_nm,convention)
    ssource=relative_daylight_spectrum(source_zenith,0.,source_zenith,turbidity,wavelength_nm,convention)
    s90=relative_daylight_spectrum(np.pi/2-source_zenith,np.pi/2,source_zenith,turbidity,wavelength_nm,convention)
    return sd,s90,ssource,spectral_factor(sd,s90,ssource)


def source_meridian_pattern(theta,phi,source_zenith,source_azimuth,separation_deg):
    """Explicit comparison convention: evaluate eq2/3 with source at +y, rotate spin2.

    This is NOT a silent change to printed eq4/5. It enforces global-plane
    rotation covariance, and still needs basis validation against camera data.
    """
    rotation=source_azimuth-np.pi/2
    return complex_pattern(theta,phi-rotation,source_zenith,np.pi/2,separation_deg)*np.exp(2j*rotation)
