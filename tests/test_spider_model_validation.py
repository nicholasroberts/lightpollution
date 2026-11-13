import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from validate_spider_sky_model import axial_residual_deg,agreement,envelope

class ValidationTests(unittest.TestCase):
    def test_axial_wrap(self):
        a=np.exp(2j*np.deg2rad(np.array([89.,-89.])))
        b=np.exp(2j*np.deg2rad(np.array([-89.,89.])))
        np.testing.assert_allclose(axial_residual_deg(a,b),[-2,2],atol=1e-10)
    def test_offset_not_hidden(self):
        b=.5*np.exp(2j*np.linspace(-1,1,40));a=b*np.exp(2j*np.deg2rad(25))
        s=agreement(a,b)
        self.assertAlmostEqual(s['raw_axial_mae_deg'],25)
        self.assertAlmostEqual(s['diagnostic_fitted_offset_deg'],25)
        self.assertLess(s['offset_removed_axial_mae_deg'],1e-10)
    def test_dilution_is_descriptive(self):
        m=np.array([0.,.2,.4,.8]);s=envelope(m*.6,m)
        self.assertEqual(s['pixels'],3)
        np.testing.assert_allclose(s['ratio_quantiles_10_50_90'],[.6]*3)
        self.assertEqual(s['fraction_above_envelope'],0)
    def test_empty(self):
        self.assertEqual(envelope(np.array([.1]),np.array([0.])),{'pixels':0})

if __name__=='__main__':unittest.main()
