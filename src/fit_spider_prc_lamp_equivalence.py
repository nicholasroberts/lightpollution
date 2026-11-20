"""Fit HPS/LED energy to reconstructed PRC; camera-green target is an RH1 proxy."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from spider_spectral_lunar import pigment_proxy


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def required_energy(target,natural_i,polarized,b,g,lamp_response,pixels):
    if not np.isfinite([target,natural_i,polarized,b,g,lamp_response,pixels]).all() or min(target,natural_i,b,g,lamp_response,pixels)<=0 or polarized<0:
        raise ValueError('Invalid target, Stokes integrals or spectral coefficients')
    natural_prc=2*np.arctanh(g*polarized/(b*natural_i))
    if target>natural_prc+1e-12:
        return dict(status='unreachable_by_unpolarized_addition',natural_prc=float(natural_prc),energy=None)
    amount=max(0.,(g*polarized/np.tanh(target/2)-b*natural_i)/(b*lamp_response*pixels))
    achieved=2*np.arctanh(g*polarized/(b*(natural_i+amount*lamp_response*pixels)))
    return dict(status='fitted',natural_prc=float(natural_prc),energy=float(amount),achieved_prc=float(achieved))


def unique_spectrum(root,filename):
    found=sorted(root.rglob(filename))
    if not found:raise ValueError('Missing '+filename+' below '+str(root))
    if len({digest(p) for p in found})!=1:raise ValueError('Different versions of '+filename+' found; pass a narrower --spectra-root')
    return found[0]


def run(a):
    if a.output.exists():raise ValueError('Use fresh output directory')
    fitpath=a.fit/'twilight_moon_fit.json';fit=json.loads(fitpath.read_text())
    targetpath=a.prc/'prc_comparison.json';target=json.loads(targetpath.read_text())
    if target['site']!=fit['site'] or target['source_fit_sha256']!=digest(fitpath):raise ValueError('PRC target and fit provenance differ')
    arrayspath=a.fit/'twilight_moon_fit_arrays.npz'
    if digest(arrayspath)!=target['source_fit_arrays_sha256']:raise ValueError('Fit arrays changed')
    if target['camera_channel']!='green':raise ValueError('Requires the agreed green-channel target')
    hps=unique_spectrum(a.spectra_root,'HPS_visible_input.csv');led=unique_spectrum(a.spectra_root,'LED_mid_visible_input.csv')
    spectra={key:np.genfromtxt(p,delimiter=',',names=True) for key,p in [('HPS',hps),('LED',led)]}
    camera=np.genfromtxt(a.camera_response,delimiter=',',names=True)
    lo=max(352.,camera['wavelength_nm'][0],*[s['wavelength_nm'][0] for s in spectra.values()])
    hi=min(700.,camera['wavelength_nm'][-1],*[s['wavelength_nm'][-1] for s in spectra.values()])
    if hi-lo<300:raise ValueError('Insufficient shared spectral band')
    wave=np.linspace(lo,hi,int(np.ceil(hi-lo))+1)
    def norm(s,column):
        values=np.interp(wave,s['wavelength_nm'],s[column])
        if not np.isfinite(values).all() or np.any(values<0) or trapezoid(values,wave)<=0:raise ValueError('Invalid energy spectrum')
        return values/trapezoid(values,wave)
    natural=norm(spectra['HPS'],'natural')
    if not np.allclose(natural,norm(spectra['LED'],'natural'),rtol=1e-6,atol=1e-10):raise ValueError('HPS/LED natural spectra differ')
    lamps={k:norm(s,'alan') for k,s in spectra.items()}
    green=np.interp(wave,camera['wavelength_nm'],camera['green_photon_response'])
    green=green/green.max()
    responses={'RH1_510':wave*pigment_proxy(wave),'D850_green_sensitivity':wave*green}
    with np.load(a.prc/'prc_support.npz',allow_pickle=False) as z:mask=np.array(z['common_central_field'])
    if mask.dtype!=np.bool_:raise ValueError('Support must be boolean')
    with np.load(arrayspath,allow_pickle=False) as z:
        if z['mixed_I'].shape!=mask.shape:raise ValueError('Support and fit dimensions differ')
        # Remove the previously fitted background completely. Only r is retained.
        natural_I=z['mixed_I']-fit['models']['mixed']['parameters'][1]
        polarization=z['mixed_Q']+1j*z['mixed_U']
    count=int(mask.sum());ni=float(natural_I[mask].sum());nz=abs(np.sum(polarization[mask]))
    if not count or ni<=0 or not np.isfinite(nz) or np.any(natural_I[mask]<=0):raise ValueError('Invalid natural Stokes support')
    target_value=target['results']['common_central_field']['measured']['maximum_absolute_prc']
    ap=target['absorbance_parallel'];at=target['absorbance_perpendicular']
    b=1-(10**(-ap)+10**(-at))/2;g=(10**(-at)-10**(-ap))/2
    results={}
    for label,response in responses.items():
        cn=float(trapezoid(natural*response,wave));coeff={k:float(trapezoid(v*response,wave)) for k,v in lamps.items()}
        fits={k:required_energy(target_value,ni*cn,nz*cn,b,g,c,count) for k,c in coeff.items()}
        ah=fits['HPS']['energy'];al=fits['LED']['energy']
        ratio=al/ah if ah is not None and ah>0 and al is not None else None
        results[label]=dict(fits=fits,natural_response_per_band_energy=cn,lamp_response_per_band_energy=coeff,
            led_to_hps_energy_ratio=ratio,led_energy_reduction_percent=None if ratio is None else 100*(1-ratio))
    report=dict(site=fit['site'],status='PRC_equivalent_spectral_scenarios_with_green_as_RH1_target_approximation',
        target_maximum_absolute_prc=target_value,target_source='new reconstructed camera-green PRC, accepted RH1 approximation',
        energy_definition='Uniform added sky radiant energy over shared band, divided by lunar zenith radiant energy over same band; not electrical power, lux or VIIRS radiance',
        spectral_band_nm=[float(lo),float(hi)],solar_to_lunar_ratio_retained=fit['models']['mixed']['parameters'][0],
        results=results,common_pixels=count,source_prc_sha256=digest(targetpath),source_support_sha256=digest(a.prc/'prc_support.npz'),
        source_fit_sha256=digest(fitpath),source_fit_arrays_sha256=digest(arrayspath),camera_response_sha256=digest(a.camera_response),
        source_spectra={k:dict(path=str(p),sha256=digest(p)) for k,p in [('HPS',hps),('LED',led)]},script_sha256=digest(__file__),
        limitations=['Sun and Moon use the same transferred natural spectral template in this first spectral comparison; the retained camera ratio then equals the assumed energy ratio.',
        'Original Sun/Moon geometry and DoLP-derived source ratio are fixed; spatial AoP was checked separately, not refitted here.',
        'Previous background is discarded; HPS/LED amounts are solved from PRC, not from DoLP.',
        'Green target approximates RH1; D850 sensitivity run repeats the forward calculation at that same numerical target, not a unique spectral reconstruction.',
        'Generic mid-CCT LED and HPS sky-transfer scenarios; no actual site lamp mix or atmospheric retrieval.',
        'ELF camera response excludes the added polarizer; published response is a sensitivity comparison.',
        'Under uniform unpolarized addition and fixed natural field, the LED/HPS energy ratio depends on spectral overlaps and is site-independent for a shared band and receptor. Site PRC sets each absolute amount.'])
    a.output.mkdir(parents=True)
    (a.output/'prc_lamp_equivalence.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    axes[0].plot(wave,pigment_proxy(wave),label='RH1 510 nm proxy');axes[0].plot(wave,green,label='ELF D850 green')
    axes[0].set(xlabel='Wavelength (nm)',ylabel='Relative photon sensitivity');axes[0].legend()
    for k,v in results.items():
        axes[1].plot([0,1],[v['fits'][lamp]['energy'] if v['fits'][lamp]['energy'] is not None else np.nan for lamp in ['HPS','LED']],'o-',label=k)
    axes[1].set_xticks([0,1]);axes[1].set_xticklabels(['HPS','LED']);axes[1].set_ylabel('Added band energy / lunar zenith band energy');axes[1].legend()
    fig.suptitle(f"{fit['site']} — same target maximum |PRC| = {target_value:.4f}\nGreen-derived target as RH1 approximation; fixed natural geometry")
    fig.savefig(a.output/'prc_lamp_equivalence.png',dpi=160);plt.close(fig)
    print(json.dumps(report['results'],indent=2));print('Saved',a.output)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fit',type=Path,required=True);p.add_argument('--prc',type=Path,required=True)
    p.add_argument('--spectra-root',type=Path,required=True)
    p.add_argument('--camera-response',type=Path,default=Path(__file__).resolve().parents[1]/'data/reference/d850_elf_response.csv')
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:run(a)
    except (OSError,ValueError,KeyError) as exc:p.error(str(exc))


if __name__=='__main__':main()
