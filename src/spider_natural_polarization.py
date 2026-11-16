"""Shared natural fractional polarization for the spider forward scenarios.

rayleigh_normalized is the investigator-selected working hypothesis, retaining
Berry phase/neutral points and the Wang horizon factor. It is not full eq17.
berry_envelope preserves the previous unit-peak phenomenological amplitude law.
Atmospheric attenuation and ALAN are applied separately by callers.
"""
import numpy as np
from wang2016_sky import berry_normalized_pattern,horizon_factor,scattering_angle
NATURAL_DOP_LAWS=('rayleigh_normalized','berry_envelope')


def natural_fractional_polarization(theta,phi,source_zenith,source_azimuth,separation_deg,law='rayleigh_normalized'):
    if law not in NATURAL_DOP_LAWS:raise ValueError('Unknown natural DoLP law')
    w=berry_normalized_pattern(theta,phi,source_zenith,source_azimuth,separation_deg)*horizon_factor(theta)
    if law=='rayleigh_normalized':
        gamma=scattering_angle(theta,phi,source_zenith,source_azimuth)
        w=w/(1+np.cos(gamma)**2)
    return w
