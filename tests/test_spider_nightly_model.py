import sys,unittest,tempfile
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from model_spider_nightly_prc import body_field,integrated_pattern,response_contrast,spectral_coefficients

class NightModelTests(unittest.TestCase):
    def test_field_and_transport(self):
        for tilt in [0,30,60]:
            rays,axis=body_field(35,tilt)
            np.testing.assert_allclose(np.linalg.norm(rays,axis=1),1,atol=1e-12)
            np.testing.assert_allclose(np.linalg.norm(axis,axis=1),1,atol=1e-12)
            np.testing.assert_allclose(np.sum(rays*axis,axis=1),0,atol=1e-12)
            self.assertGreater(rays[:,2].min(),0)
    def test_world_rotation_invariance(self):
        values=[]
        for h,az in [(20,40),(85,105)]:
            rays,axis=body_field(h,30,4096)
            values.append(integrated_pattern(rays,axis,25,az,24))
        self.assertAlmostEqual(*values,places=10)
    def test_alan_and_atmosphere_reduce_contrast(self):
        coefficients,_=spectral_coefficients()
        full=abs(response_contrast(.5,1,0,coefficients))
        self.assertLess(abs(response_contrast(.5,1,1,coefficients)),full)
        self.assertLess(abs(response_contrast(.5,.8,0,coefficients)),full)
        self.assertEqual(response_contrast(.5,0,3,coefficients),0)
    def test_spectral_normalization_and_weighting(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'spectra.csv'
            path.write_text('wavelength_nm,natural,alan,sensitivity\n400,1,3,1\n500,2,2,1\n600,3,1,1\n')
            actual,_=spectral_coefficients(path);expected,_=spectral_coefficients()
            np.testing.assert_allclose(actual,expected)
