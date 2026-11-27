"""Produce Figure 2A–C using saved Jones lunar spectra and locked Figure 1 fits."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_paper1_fig02a import spectral_inputs
from paper1_jones_spectral_input import load_jones_moon
from paper1_figure_archive import archive_figure,digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--jones-root',type=Path,default=Path('data/processed/paper1_jones_moon_spectra_01'))
    p.add_argument('--spectra-root',type=Path,default=Path('data/processed/paper1_spectral_inputs_02'))
    p.add_argument('--fit-root',type=Path,default=Path('data/processed/paper1_fig01c_joint_separation_01'))
    p.add_argument('--figure-root',type=Path,default=Path('figures/paper1'))
    p.add_argument('--figure-version',default='v04')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if Path(a.figure_version).name!=a.figure_version or a.figure_version in ('.','..'):
        p.error('Figure version must be a directory name')
    destinations={panel:a.figure_root/('fig02'+panel)/a.figure_version for panel in ('a','b','c')}
    if a.output.exists() or any(path.exists() for path in destinations.values()):
        p.error('Use fresh output and figure-version destinations')
    wave,spectra,source_report=spectral_inputs(a.spectra_root)
    moons={};records={}
    for site in ('Cottage','museum'):
        jw,jm,metadata=load_jones_moon(a.jones_root,site,a.fit_root/site/'fit/solar_lunar_fit.json')
        if not np.array_equal(jw,wave):p.error('Jones and lamp spectral grids differ')
        moons[site]=jm;records[site]=metadata
    a.output.mkdir(parents=True)
    panel=a.output/'fig02a';panel.mkdir()
    plt.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(7.5,4.3),layout='constrained')
    series=[
        ('Jones moonlight: Cottage',moons['Cottage'],'#0072B2','-'),
        ('Jones moonlight: museum',moons['museum'],'#56B4E9','--'),
        ('Retained solar spectral proxy',spectra['Moon'],'#777777',':'),
        ('HPS',spectra['HPS'],'#E69F00','-'),
        ('LED',spectra['LED'],'#009E73','-')]
    for label,values,color,style in series:
        ax.plot(wave,values,color=color,ls=style,lw=1.5,label=label)
    ax.set(xlabel='Wavelength (nm)',ylabel='Normalized spectral energy (nm⁻¹)',
        xlim=(400,700),ylim=(0,None),title='Natural and ALAN spectral inputs\n400–700 nm; equal integrated radiant energy')
    ax.legend(frameon=False,fontsize=8);ax.spines[['top','right']].set_visible(False)
    fig.text(.015,.99,'A',ha='left',va='top',fontsize=18,fontweight='bold')
    for ext in ('png','pdf','svg'):fig.savefig(panel/('fig02a_spectra.'+ext),dpi=300)
    plt.close(fig)
    np.savetxt(panel/'spectra_equal_band_energy.csv',
        np.column_stack([wave,moons['Cottage'],moons['museum'],spectra['Moon'],spectra['HPS'],spectra['LED']]),
        delimiter=',',comments='',header='wavelength_nm,Jones_Moon_Cottage,Jones_Moon_museum,solar_proxy,HPS,LED')
    report=dict(subfigure='2A',spectral_band_nm=[400,700],
        normalization='Unit integrated energy for each curve; not photon, photopic or electrical units',
        jones_inputs=records,retained_solar_and_lamp_inputs=source_report,
        script_sha256=digest(__file__),solar_spectrum='Retained prior transferred natural spectrum as approximate solar proxy; not an independent twilight calculation',
        atmosphere='Jones lunar spectra use fixed ESO reference atmosphere; no additional transfer applied',
        field_approximation='Zenith lunar spectral shape used throughout 60-degree field; Figure 1 camera-green mixture preserved')
    (panel/'spectral_inputs.json').write_text(json.dumps(report,indent=2)+'\n')
    archive_figure(panel,destinations['a'])
    src=Path(__file__).resolve().parent;summary={}
    for site,letter in [('Cottage','b'),('museum','c')]:
        output=a.output/site/'fit'
        subprocess.run([sys.executable,str(src/('fit_paper1_fig02'+letter+'.py')),
            '--fit',str(a.fit_root/site/'fit'),'--spectra-root',str(a.spectra_root),
            '--jones-root',str(a.jones_root),'--output',str(output),
            '--figure-copy',str(destinations[letter])],check=True)
        fit=json.loads((output/'lamp_DoLP_fit.json').read_text())
        summary[site]=dict(solar_fraction_camera_green=fit['solar_fraction'],
            natural_band_weights=fit['natural_band_weights'],
            measured_mean_DoLP=fit['target_mean_camera_green_dolp'],
            fits={lamp:{key:values[key] for key in ('amount','achieved_mean_dolp','spatial_DoLP_rms_residual','mean_DoLP_residual','optimizer_success')}
                for lamp,values in fit['fits'].items()},
            fit_report_path=str(output/'lamp_DoLP_fit.json'),
            fit_report_sha256=digest(output/'lamp_DoLP_fit.json'))
    (a.output/'fig02_jones_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    print('Archived Figure 2A–C:',', '.join(str(v) for v in destinations.values()))

if __name__=='__main__':main()
