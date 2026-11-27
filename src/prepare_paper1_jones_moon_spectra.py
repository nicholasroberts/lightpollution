"""Jones/ESO scattered-moonlight reference spectra at each capture's zenith.

A simple spectral comparison only: fixed Paranal reference atmosphere, no fits,
no change to Figure 1 or the production ALAN inputs.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from astropy.io import fits
from astropy.coordinates import get_body,get_sun
from astropy.time import Time
from astropy.utils import iers
from paper1_figure_archive import archive_figure

H=6.62607015e-34
C=299792458.
ARCSEC2_PER_SR=(180*3600/np.pi)**2


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def integral(y,x):
    return np.trapezoid(y,x) if hasattr(np,'trapezoid') else np.trapz(y,x)


def energy_from_photons(wavelength_nm,photons):
    """SkyCalc ph/s/m2/micron/arcsec2 -> W/m2/nm/sr."""
    w=np.asarray(wavelength_nm,dtype=float)
    return np.asarray(photons,dtype=float)*H*C/(w*1e-9)/1000.*ARCSEC2_PER_SR


def unit_energy(wave,energy):
    area=integral(energy,wave)
    if not np.isfinite(area) or area<=0:raise ValueError('Nonpositive moonlight energy')
    return energy/area


def read_moon_component(path):
    with fits.open(path) as hdul:
        table=hdul[1]
        names={k.upper():k for k in table.columns.names}
        if not {'LAM','FLUX_SML'}<=set(names):
            raise ValueError('SkyCalc FITS requires LAM and FLUX_SML, not total FLUX')
        column=table.columns[names['LAM']]
        unit=(column.unit or '').strip().lower()
        if unit not in ('nm','nanometer','nanometre'):
            raise ValueError('Expected documented SkyCalc 1.6 wavelength unit nm; found '+repr(unit))
        wave=np.array(table.data[names['LAM']],dtype=float)
        photons=np.array(table.data[names['FLUX_SML']],dtype=float)
        flux_unit=table.columns[names['FLUX_SML']].unit
        if flux_unit!='ph/s/m2/micron/arcsec2':
            raise ValueError('Unexpected SkyCalc photon radiance unit: '+repr(flux_unit))
    if wave.ndim!=1 or wave.size<2 or not np.isfinite(wave).all() or np.any(np.diff(wave)<=0):
        raise ValueError('Invalid wavelength grid')
    if not np.isfinite(photons).all() or np.any(photons<0):
        raise ValueError('Invalid scattered-moonlight photon radiance')
    if wave[0]>400.+1e-6 or wave[-1]<700.-1e-6:
        raise ValueError('Moonlight output must cover 400–700 nm without extrapolation')
    target=np.arange(400.,701.)
    photon_interp=np.interp(target,wave,photons)
    energy=energy_from_photons(target,photon_interp)
    return target,photon_interp,energy,dict(wavelength_unit=unit,photon_radiance_unit=flux_unit)


def parameters(anchor):
    altitude=float(anchor['moon_elevation_deg'])
    if not 0<altitude<90:raise ValueError('Moon must be above horizon for this simple comparison')
    iers.conf.auto_download=False
    moment=Time(anchor['utc'].replace('Z',''),scale='utc')
    moon=get_body('moon',moment);sun=get_sun(moment)
    # Geocentric elongation is the documented SkyCalc moon_sun_sep input.
    elongation=float(moon.separation(sun).deg)
    distance=float(moon.distance.to_value('km')/384400.)
    return dict(airmass=1.,moon_alt=altitude,moon_target_sep=90.-altitude,
        moon_sun_sep=elongation,moon_earth_dist=distance,
        incl_moon='Y',incl_starlight='N',incl_zodiacal='N',
        incl_loweratm='N',incl_upperatm='N',incl_airglow='N',incl_therm='N',
        pwv_mode='pwv',pwv=3.5,observatory='paranal',
        vacair='vac',wmin=400.,wmax=700.,wgrid_mode='fixed_wavelength_step',
        wdelta=1.,lsf_type='none',temp_flag=0)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fit-root',type=Path,default=Path('data/processed/paper1_fig01c_joint_separation_01'))
    p.add_argument('--sites',nargs='+',default=['Cottage','museum'])
    p.add_argument('--spectra-root',type=Path,default=Path('data/processed/paper1_spectral_inputs_02'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--figure-copy',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists() or a.figure_copy.exists():p.error('Use fresh output and figure-copy paths')
    # Check source coverage before creating output or making remote requests.
    source=a.spectra_root/'moon_source_352_750.csv'
    lunar=np.genfromtxt(source,delimiter=',',names=True)
    sw,se=lunar['wavelength_nm'],lunar['relative_energy']
    if not np.isfinite(sw).all() or np.any(np.diff(sw)<=0) or not np.isfinite(se).all() or np.any(se<0) or sw[0]>400 or sw[-1]<700:
        p.error('Audited lunar source must cover 400–700 nm without extrapolation')
    # Preparation stores the source audit in inherited spectral provenance.
    inherited=json.loads((a.spectra_root/'spectral_provenance.json').read_text())
    if inherited.get('moon_audit',{}).get('status')!='PASSED_COARSE_SPECTRAL_SHAPE_CHECK_NOT_SITE_VALIDATION':
        p.error('Requires the wider audited Moon preparation, not the old NIST table')
    paths={site:a.fit_root/site/'fit/solar_lunar_fit.json' for site in a.sites}
    geometry={}
    for site,path in paths.items():
        report=json.loads(path.read_text())
        if report['site']!=site:p.error('Fit site mismatch')
        geometry[site]=report['astronomical_anchor']
    version=importlib.metadata.version('skycalc_cli')
    if tuple(map(int,version.split('.')[:2]))<(1,6):p.error('Install skycalc_cli>=1.6')
    a.output.mkdir(parents=True)
    results={};curves={}
    for site in a.sites:
        folder=a.output/site;folder.mkdir()
        inputs=parameters(geometry[site]);input_file=folder/'skycalc_input.json'
        input_file.write_text(json.dumps(inputs,indent=2)+'\n')
        output_file=folder/'skycalc_moon.fits'
        print('Requesting Jones/ESO zenith spectrum: '+site,flush=True)
        # Official CLI/API: no handwritten replacement for Jones scattering.
        subprocess.run([sys.executable,'-c','from skycalc_cli.skycalc_cli import main; main()',
            '-i',str(input_file),'-o',str(output_file)],check=True,timeout=180)
        if not output_file.is_file():raise RuntimeError('SkyCalc did not produce its FITS output')
        wave,photons,energy,units=read_moon_component(output_file)
        source_shape=unit_energy(wave,np.interp(wave,sw,se))
        old_shape=unit_energy(wave,np.interp(wave,sw,se)*(550./wave)**4)
        shape=unit_energy(wave,energy)
        np.savetxt(folder/'moon_spectra_400_700.csv',
            np.column_stack([wave,photons,energy,shape,source_shape,old_shape]),
            delimiter=',',comments='',
            header='wavelength_nm,jones_photons_per_s_m2_micron_arcsec2,jones_energy_W_m2_nm_sr,jones_unit_band_energy,lunar_source_unit_band_energy,old_rayleigh_unit_band_energy')
        curves[site]=(wave,source_shape,old_shape,shape)
        results[site]=dict(capture_anchor=geometry[site],source_fit_path=str(paths[site]),
            source_fit_sha256=digest(paths[site]),skycalc_parameters=inputs,
            skycalc_input_sha256=digest(input_file),skycalc_fits_sha256=digest(output_file),
            fits_units=units,jones_peak_sample_nm=float(wave[np.argmax(shape)]),
            old_peak_sample_nm=float(wave[np.argmax(old_shape)]),
            jones_band_energy_W_m2_sr=float(integral(energy,wave)))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,len(a.sites),figsize=(6*len(a.sites),4),squeeze=False,layout='constrained')
    for ax,site in zip(axes.ravel(),a.sites):
        w,src,old,new=curves[site]
        ax.plot(w,src,color='#777777',label='Lunar source (before atmosphere)')
        ax.plot(w,old,color='#E69F00',ls='--',label='Previous source × wavelength⁻⁴')
        ax.plot(w,new,color='#0072B2',label='Jones/ESO scattered moonlight')
        ax.set(xlabel='Wavelength (nm)',ylabel='Normalized spectral energy (nm⁻¹)',
            xlim=(400,700),ylim=(0,None),title=site+'; Moon '+f"{geometry[site]['moon_elevation_deg']:.1f}°; view: zenith")
        ax.legend(frameon=False,fontsize=8)
    fig.suptitle('Scattered moonlight comparison — fixed ESO reference atmosphere\n400–700 nm; equal band energy; no spectral fit')
    panel=a.output/'panel';panel.mkdir()
    for ext in ('png','pdf','svg'):fig.savefig(panel/('jones_moon_spectral_comparison.'+ext),dpi=220)
    plt.close(fig)
    report=dict(status='JONES_ESO_REFERENCE_ATMOSPHERE_COMPARISON_NOT_LOWLAND_SITE_CALIBRATION',
        model='Official ESO SkyCalc implementation of Jones et al. 2013',
        model_doi='10.1051/0004-6361/201322433',skycalc_cli_version=version,
        comparison_source_moon_path=str(source),comparison_source_moon_sha256=digest(source),
        source_provenance=inherited,script_sha256=digest(__file__),sites=results,
        spectral_band_nm=[400,700],view='Zenith only; not a field-average spectrum',
        atmosphere='Fixed stock Paranal reference, PWV 3.5 mm; no fitted atmospheric parameters',
        units='FLUX_SML photon radiance converted once via hc/lambda to energy, then normalized',
        figures_and_ALAN_fits_modified=False,
        limitations=['Reference atmosphere is not calibrated for lowland Cottage/museum.',
            'No full-field spectral map or polarization radiative transfer is claimed.',
            'No multiplication by an additional wavelength^-4 or TRANS: SkyCalc output already includes propagation.',
            'No solar/twilight spectrum supplied: do not use Moon output for museum solar contribution.',
            'Direct source and scattered sky curves represent different quantities; source shown for context.',
            'Jones/ESO uses its own solar/lunar source prescription; the supplied lunar CSV is only the context and old-transfer reference.',
            'A bluer or redder peak is an outcome, not a fitted requirement.'])
    for root in (a.output,panel):(root/'jones_moon_spectra.json').write_text(json.dumps(report,indent=2)+'\n')
    archive_figure(panel,a.figure_copy)
    print('Saved spectral comparison:',a.figure_copy)

if __name__=='__main__':main()
