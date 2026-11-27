"""Read a saved official Jones/ESO lunar spectrum with provenance checks."""
import json
from pathlib import Path
import numpy as np
from prepare_paper1_jones_moon_spectra import read_moon_component,unit_energy
from paper1_figure_archive import digest


def load_jones_moon(root,site,expected_fit=None):
    root=Path(root);rp=root/'jones_moon_spectra.json'
    report=json.loads(rp.read_text());entry=report['sites'][site]
    if report.get('model_doi')!='10.1051/0004-6361/201322433':
        raise ValueError('Jones/ESO report required')
    if expected_fit is not None and entry['source_fit_sha256']!=digest(expected_fit):
        raise ValueError('Jones capture and final Figure 1 fit differ: '+site)
    fits_path=root/site/'skycalc_moon.fits'
    if digest(fits_path)!=entry['skycalc_fits_sha256']:
        raise ValueError('Jones FITS differs from recorded output: '+site)
    wave,photons,energy,units=read_moon_component(fits_path)
    return wave,unit_energy(wave,energy),dict(
        report_path=str(rp),report_sha256=digest(rp),
        fits_path=str(fits_path),fits_sha256=digest(fits_path),
        capture_anchor=entry['capture_anchor'],parameters=entry['skycalc_parameters'],
        atmosphere=report['atmosphere'],view=report['view'],
        model_doi=report['model_doi'],conversion='Photon radiance converted once to energy; no further scattering or extinction',
        solar_treatment='Jones replaces lunar component only; existing transferred natural spectrum retained as approximate solar proxy')


def band_weights_from_green(fraction,moon_overlap,solar_overlap):
    """Preserve fitted camera-green component weights after changing spectra."""
    if not 0<=fraction<=1 or min(moon_overlap,solar_overlap)<=0:
        raise ValueError('Invalid green fraction or spectral overlap')
    weights=np.array([(1-fraction)/moon_overlap,fraction/solar_overlap])
    return weights/weights.sum()
