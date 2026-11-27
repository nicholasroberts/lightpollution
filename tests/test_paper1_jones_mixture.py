"""Preserve the approved green pattern when lunar/solar spectra differ."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from paper1_jones_spectral_input import band_weights_from_green

class JonesMixtureTests(unittest.TestCase):
    def test_camera_green_stokes_shape_is_preserved(self):
        f=.743078968425043
        km,ks=180.,245.
        wm,ws=band_weights_from_green(f,km,ks)
        self.assertAlmostEqual(wm+ws,1.)
        self.assertAlmostEqual(ws*ks/(wm*km+ws*ks),f)
        moon=np.array([[1.,.2,-.1],[1.,-.3,.2]])
        solar=np.array([[1.,-.1,.4],[1.,.25,-.2]])
        old=(1-f)*moon+f*solar
        new=wm*km*moon+ws*ks*solar
        np.testing.assert_allclose(new/new[:,0,None],old/old[:,0,None],rtol=1e-13)

    def test_identical_spectra_recover_original_band_weights(self):
        np.testing.assert_allclose(band_weights_from_green(.2,200.,200.),[.8,.2])

    def test_lunar_and_solar_endpoints(self):
        np.testing.assert_array_equal(band_weights_from_green(0.,100.,200.),[1.,0.])
        np.testing.assert_array_equal(band_weights_from_green(1.,100.,200.),[0.,1.])
        with self.assertRaises(ValueError):band_weights_from_green(.2,0.,200.)

if __name__=='__main__':unittest.main()
