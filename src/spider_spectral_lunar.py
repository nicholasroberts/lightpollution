"""Declared pigment proxy and KS1991 V-band relative-brightness pilot.

KS1991 doi10.1086/132921 eqs3,8,9,15,21; separation >=10deg.
Stavenga alpha template: x=log10(lambda/lambda_max), exp(-380*x^2*(1+6.09*x+3/8*(6.09*x)^2)).
Default is one green-sensitive alpha band at510nm. Optional legacy UV Gaussian
is inactive by default; no molecular RH1 identity is inferred.
"""
import numpy as np


def pigment_proxy(wavelength_nm,peak_nm=510.,uv_peak_nm=365.,uv_amplitude=0.,uv_width_nm=40.):
    w=np.asarray(wavelength_nm,dtype=float)
    if np.any(~np.isfinite(w)) or np.any(w<=0) or not np.isfinite([peak_nm,uv_peak_nm,uv_amplitude,uv_width_nm]).all() or min(peak_nm,uv_peak_nm,uv_width_nm)<=0 or uv_amplitude<0:raise ValueError('Invalid pigment parameters')
    x=np.log10(w/peak_nm);bx=6.09*x
    alpha=np.exp(-380*x*x*(1+bx+3/8*bx*bx))
    uv=uv_amplitude*np.exp(-((w-uv_peak_nm)/uv_width_nm)**2)
    # Common visible-peak normalization, not normalization on caller's grid.
    reference=1+uv_amplitude*np.exp(-((peak_nm-uv_peak_nm)/uv_width_nm)**2)
    return (alpha+uv)/reference


def lunar_phase_distance(utc_times):
    from astropy.coordinates import get_body,get_sun,solar_system_ephemeris
    from astropy.time import Time
    from astropy.utils import iers
    iers.conf.auto_download=False
    from datetime import datetime,timezone
    parsed=[datetime.fromisoformat(v.replace('Z','+00:00')) if isinstance(v,str) else v for v in utc_times]
    parsed=[v.replace(tzinfo=timezone.utc) if v.tzinfo is None else v.astimezone(timezone.utc) for v in parsed]
    t=Time(parsed)
    with solar_system_ephemeris.set('builtin'):
        moon=get_body('moon',t).cartesian.xyz.to_value('km').T
        sun=get_sun(t).cartesian.xyz.to_value('km').T
    tosun=sun-moon;toearth=-moon
    cosine=np.sum(tosun*toearth,axis=-1)/(np.linalg.norm(tosun,axis=-1)*np.linalg.norm(toearth,axis=-1))
    return np.rad2deg(np.arccos(np.clip(cosine,-1,1))),np.linalg.norm(moon,axis=-1)


def ks91_brightness(view_zenith_rad,moon_elevation_deg,separation_rad,phase_angle_deg,extinction=.2,distance_km=384400.):
    """Moon-only V-band brightness in nL, plus explicit inverse-square distance factor.

    Not a spectral model; exclusion below horizon and <10deg is deliberate.
    Extinction controls I, independent of phenomenological Q/U depolarization.
    """
    theta=np.asarray(view_zenith_rad,dtype=float);rho=np.asarray(separation_rad,dtype=float)
    if not np.isfinite([moon_elevation_deg,phase_angle_deg,extinction,distance_km]).all() or not 0<moon_elevation_deg<90 or not 0<=phase_angle_deg<=180 or extinction<=0 or distance_km<=0:raise ValueError('Invalid lunar geometry/extinction')
    if np.any(~np.isfinite(theta)) or np.any((theta<0)|(theta>np.pi/2)) or np.any(~np.isfinite(rho)) or np.any((rho<np.deg2rad(10))|(rho>np.pi)):raise ValueError('KS91 pilot requires upper sky and separation >=10deg')
    xm=(1-.96*np.cos(np.deg2rad(moon_elevation_deg))**2)**(-.5)
    xv=(1-.96*np.sin(theta)**2)**(-.5)
    a=abs(phase_angle_deg)
    source=10**(-.4*(3.84+.026*a+4e-9*a**4))*(384400./distance_km)**2
    scatter=10**5.36*(1.06+np.cos(rho)**2)+10**(6.15-np.rad2deg(rho)/40)
    return scatter*source*10**(-.4*extinction*xm)*(1-10**(-.4*extinction*xv))
