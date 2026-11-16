"""Isolated normalization diagnostic: not a replacement for the Wang/KS91 pilot.

No measured-site fitting, production settings or biological predictions changed.
"""
import argparse,json
from pathlib import Path
import numpy as np
from wang2016_sky import berry_normalized_pattern


def ideal_rayleigh(gamma):
    g=np.asarray(gamma,dtype=float)
    return np.sin(g)**2/(1+np.cos(g)**2)


def mixed_dop(natural_intensity,natural_dop,alan_intensity):
    i,p,a=np.broadcast_arrays(np.asarray(natural_intensity,dtype=float),np.asarray(natural_dop,dtype=float),np.asarray(alan_intensity,dtype=float))
    if np.any(~np.isfinite(i+p+a)) or np.any(i<=0) or np.any((p<0)|(p>1)) or np.any(a<0):raise ValueError('Invalid physical Stokes inputs')
    return i*p/(i+a)


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    gamma=np.linspace(0,np.pi,18001);q=np.sin(gamma)**2;i=1+np.cos(gamma)**2;dop=q/i
    # Horizon source/equatorial directions isolate scattering angle. No horizon
    # attenuation is applied: this is an algebraic limit, not an actual sky map.
    kernel=abs(berry_normalized_pattern(np.full_like(gamma,np.pi/2),np.pi/2+gamma,np.pi/2,np.pi/2,.001))
    widths=dict(polarized_intensity_fwhm_deg=90.,rayleigh_dop_fwhm_deg=float(2*np.rad2deg(np.arcsin(1/np.sqrt(3)))),
        rayleigh_plus_absolute_alan_fwhm_deg={str(a):float(2*np.rad2deg(np.arcsin(np.sqrt((1+a)/(3+2*a)))) ) for a in [0.,.5,1.,3.]})
    a.output.mkdir(parents=True)
    fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    x=np.rad2deg(gamma)
    axes[0].plot(x,kernel,label='Berry nearly coincident neutral points: |W|')
    axes[0].plot(x,dop,label='Rayleigh DoLP: sin²γ/(1+cos²γ)');axes[0].axhline(.5,color='gray',ls='--')
    axes[0].set(xlabel='Scattering angle (deg)',ylabel='Unit-peak amplitude',title='Polarized intensity versus DoLP');axes[0].legend(fontsize=8)
    for addition in [0.,.5,1.,3.]:
        values=mixed_dop(i,dop,addition)
        axes[1].plot(x,values/np.max(values),label=f'Uniform absolute ALAN={addition:g}')
    axes[1].set(xlabel='Scattering angle (deg)',ylabel='Peak-normalized DoLP',title='Consistent ideal Stokes + additive ALAN');axes[1].legend(fontsize=8)
    fig.suptitle('Normalization diagnostic only: no horizon factor, split=.001deg, no site fits\nRayleigh ideal atmosphere; not the Wang spectral model or a new animal prediction')
    fig.savefig(a.output/'intensity_normalization_audit.png',dpi=180);plt.close(fig)
    np.savetxt(a.output/'ideal_profiles.csv',np.column_stack([x,kernel,q,i,dop]),delimiter=',',header='scattering_angle_deg,berry_magnitude,ideal_polarized_intensity,ideal_total_intensity,ideal_dop',comments='')
    (a.output/'audit.json').write_text(json.dumps(dict(status='ALGEBRAIC_DIAGNOSTIC_NOT_PRODUCTION_MODEL',widths=widths,
        berry_limit_max_abs_error_vs_sin_squared=float(np.max(abs(kernel-q))),
        production_model_changed=False,sites_fitted=False,
        stokes_contract='Q+iU = I_natural*p_natural*exp(2i*psi); I_total=I_natural+I_ALAN; Q,U unchanged by unpolarized ALAN',
        interpretation=['Wang treats geometric magnitude as phenomenological DoLP; this diagnostic does not establish an omitted denominator in its implementation.',
        'Do not divide an already normalized DoLP by a second arbitrary brightness model.',
        'The right panel uses Q=sin²gamma and I=1+cos²gamma consistently. It does not combine Berry dimensional polarized intensity with an independently scaled KS91 total intensity.',
        'Constant eta or constant ALAN/natural ratio preserves normalized width. Constant absolute ALAN may change width; direction/sign depends on the natural intensity field.']),indent=2)+'\n')
    print(a.output/'intensity_normalization_audit.png')

if __name__=='__main__':main()
