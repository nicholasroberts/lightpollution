"""Numerical checks for the official Jones/ESO spectrum adapter."""
import sys
from pathlib import Path
import tempfile
import unittest
import numpy as np
from astropy.io import fits

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from prepare_paper1_jones_moon_spectra import (
    H,C,ARCSEC2_PER_SR,energy_from_photons,unit_energy,integral,
    read_moon_component,parameters)

class JonesSpectrumTests(unittest.TestCase):
    def test_photon_energy_conversion_once(self):
        w=np.array([400.,700.])
        e=energy_from_photons(w,np.ones(2))
        np.testing.assert_allclose(e,H*C/(w*1e-9)*ARCSEC2_PER_SR/1000.)
        self.assertAlmostEqual(e[0]/e[1],700./400.)
        self.assertAlmostEqual(integral(unit_energy(w,e),w),1.)

    def test_extract_moon_not_total_and_reject_wrong_units(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'sky.fits';w=np.arange(400.,701.)
            cols=[fits.Column(name='LAM',format='D',unit='nm',array=w),
                  fits.Column(name='FLUX_SML',format='D',unit='ph/s/m2/micron/arcsec2',array=np.ones_like(w)),
                  fits.Column(name='FLUX',format='D',array=100*np.ones_like(w))]
            fits.HDUList([fits.PrimaryHDU(),fits.BinTableHDU.from_columns(cols)]).writeto(path)
            wave,photons,energy,units=read_moon_component(path)
            np.testing.assert_array_equal(photons,np.ones_like(w))
            np.testing.assert_allclose(energy,energy_from_photons(w,photons))
            cols[0]=fits.Column(name='LAM',format='D',unit='micron',array=w/1000)
            fits.HDUList([fits.PrimaryHDU(),fits.BinTableHDU.from_columns(cols)]).writeto(path,overwrite=True)
            with self.assertRaisesRegex(ValueError,'unit nm'):
                read_moon_component(path)

    def test_zenith_capture_geometry_and_fixed_settings(self):
        p=parameters(dict(utc='2018-08-21T20:23:49Z',moon_elevation_deg=12.31370836118817))
        self.assertAlmostEqual(p['moon_target_sep']+p['moon_alt'],90.)
        self.assertEqual(p['airmass'],1.)
        self.assertEqual([p['wmin'],p['wmax']],[400.,700.])
        self.assertEqual(p['observatory'],'paranal')
        self.assertEqual(p['incl_airglow'],'N')
        self.assertTrue(0<p['moon_sun_sep']<180)
        self.assertTrue(.91<p['moon_earth_dist']<1.08)

if __name__=='__main__':unittest.main()
