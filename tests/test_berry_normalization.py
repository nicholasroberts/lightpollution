import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wang2016_sky import berry_normalized_pattern,source_meridian_pattern

class BerryNormalizationTests(unittest.TestCase):
    def test_prefactor_and_sign(self):
        th=np.array([.1,.4,1.]);ph=np.array([.2,.8,1.2])
        np.testing.assert_allclose(berry_normalized_pattern(th,ph,.7,.9,24),
                                 -4*source_meridian_pattern(th,ph,.7,.9,24))
    def test_analytic_unity_direction(self):
        # +/-x horizon lies 90deg from both meridional neutral points.
        for zs in [.3,.7,1.3]:
            for separation in [12,24,36]:
                self.assertAlmostEqual(abs(berry_normalized_pattern(np.pi/2,0,zs,np.pi/2,separation)),1.)
    def test_coalesced_neutrals_limit_is_sine_squared(self):
        # Berry intensity limit is sin^2(gamma), NOT Rayleigh normalized DoLP.
        th=np.linspace(.01,1.5,90);ph=np.linspace(-3,3,90);zs=.7;az=.9
        cosgamma=np.cos(th)*np.cos(zs)+np.sin(th)*np.sin(zs)*np.cos(ph-az)
        np.testing.assert_allclose(abs(berry_normalized_pattern(th,ph,zs,az,1e-5)),1-cosgamma**2,atol=1e-12)
