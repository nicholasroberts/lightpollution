"""Figure 2A: archived Moon, HPS and LED spectral inputs at equal band energy."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from paper1_figure_archive import archive_figure,digest


def integral(values,wave):
    return np.trapezoid(values,wave) if hasattr(np,'trapezoid') else np.trapz(values,wave)


def load_spectrum(path):
    s=np.genfromtxt(path,delimiter=',',names=True)
    if s.dtype.names is None or not {'wavelength_nm','natural','alan'}<=set(s.dtype.names):
        raise ValueError('Required CSV columns: wavelength_nm,natural,alan: '+str(path))
    if s.ndim!=1 or len(s)<2 or not np.isfinite(s['wavelength_nm']).all() or np.any(np.diff(s['wavelength_nm'])<=0):
        raise ValueError('Wavelengths must be finite and strictly increasing')
    for key in ('natural','alan'):
        if not np.isfinite(s[key]).all() or np.any(s[key]<0):raise ValueError('Invalid energy spectrum: '+key)
    return s


def spectral_inputs(root):
    paths={k:root/n for k,n in [('HPS','HPS_visible_input.csv'),('LED','LED_mid_visible_input.csv')]}
    data={k:load_spectrum(p) for k,p in paths.items()}
    provenance=root/'spectral_provenance.json'
    metadata=json.loads(provenance.read_text())
    # Fixed publication band: fail rather than silently narrowing old inputs.
    lo,hi=400.,700.
    for label,values in data.items():
        if values['wavelength_nm'][0]>lo or values['wavelength_nm'][-1]<hi:
            raise ValueError(label+' input must cover the full 400–700 nm band without extrapolation; regenerate wider inputs with prepare_spider_spectral_pilot.py')
    wave=np.unique(np.concatenate([np.linspace(lo,hi,int(np.ceil(hi-lo))+1)]+
        [v['wavelength_nm'][(v['wavelength_nm']>=lo)&(v['wavelength_nm']<=hi)] for v in data.values()]))
    def norm(s,key):
        v=np.interp(wave,s['wavelength_nm'],s[key]);area=integral(v,wave)
        if area<=0:raise ValueError('Spectrum has zero band energy')
        return v/area
    moon=norm(data['HPS'],'natural')
    if not np.allclose(moon,norm(data['LED'],'natural'),rtol=1e-6,atol=1e-10):
        raise ValueError('Natural spectral templates differ between HPS and LED input files')
    spectra={'Moon':moon,'HPS':norm(data['HPS'],'alan'),'LED':norm(data['LED'],'alan')}
    report=dict(spectral_band_nm=[float(lo),float(hi)],normalization='Unit integrated radiant energy within the shared band; normalized energy per nm, not photons or human luminance',
        source_files={k:dict(path=str(p),sha256=digest(p),available_band_nm=[float(data[k]['wavelength_nm'][0]),float(data[k]['wavelength_nm'][-1])]) for k,p in paths.items()},
        source_provenance_path=str(provenance),source_provenance_sha256=digest(provenance),
        inherited_spectral_provenance=metadata,
        natural_template='Existing transferred natural/Moon template shared by the HPS/LED inputs; not a new site-specific measured lunar spectrum',
        sky_transfer='Existing input transfer retained; no additional scattering/extinction applied by this plot',
        response_weighting_applied=False,fitting_performed=False,
        band_integrals={k:float(integral(v,wave)) for k,v in spectra.items()})
    return wave,spectra,report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spectra-root',type=Path,default=Path('data/processed/paper1_spectral_inputs_02'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--figure-copy',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists() or a.figure_copy.exists():p.error('Use fresh output and figure-copy directories')
    try:wave,spectra,report=spectral_inputs(a.spectra_root)
    except (OSError,ValueError) as exc:p.error(str(exc))
    print('Shared spectral band: '+f"{wave[0]:g}–{wave[-1]:g} nm"+'; each spectrum normalized over this interval',flush=True)
    plt.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
    for label,color in [('Moon','#0072B2'),('HPS','#E69F00'),('LED','#009E73')]:
        ax.plot(wave,spectra[label],label=label,color=color,lw=1.6)
    ax.set(xlabel='Wavelength (nm)',ylabel='Normalized spectral energy (nm⁻¹)',
        xlim=(wave[0],wave[-1]),ylim=(0,None),title='Moon, HPS and LED spectra\nEqual integrated radiant energy')
    ax.legend(frameon=False);ax.spines[['top','right']].set_visible(False)
    fig.text(.015,.99,'A',ha='left',va='top',fontsize=18,fontweight='bold')
    a.output.mkdir(parents=True)
    for ext in ('png','pdf','svg'):fig.savefig(a.output/('fig02a_spectra.'+ext),dpi=300)
    plt.close(fig)
    np.savetxt(a.output/'spectra_equal_band_energy.csv',np.column_stack([wave,*spectra.values()]),
        delimiter=',',header='wavelength_nm,Moon,HPS,LED',comments='')
    report.update(subfigure='2A',script_sha256=digest(__file__))
    (a.output/'spectral_inputs.json').write_text(json.dumps(report,indent=2)+'\n')
    archive_figure(a.output,a.figure_copy)
    print('Saved',a.figure_copy)


if __name__=='__main__':main()
