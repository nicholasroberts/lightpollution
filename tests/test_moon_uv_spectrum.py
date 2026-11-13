import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from prepare_moon_uv_spectrum import normalize,bin_means
from plot_spider_spectral_pilot import fov_statistics
from prepare_site_alan_scenarios import overlap_conversion

class MoonUVTests(unittest.TestCase):
    def test_shape_normalization_removes_absolute_calibration(self):
        w=np.arange(352.,751.);v=1+w/1000
        np.testing.assert_allclose(normalize(w,v),normalize(w,17*v))
        np.testing.assert_allclose(bin_means(w,w,np.array([400,410,420])),[405,415])
    def test_fov_statistics_excludes_outside_and_invalid(self):
        a=np.full((7,7),99.);a[2:5,2:5]=np.array([[.1,.2,.3],[.4,np.nan,.6],[.7,.8,.9]])
        stats=fov_statistics(a,{'center':3,'fov_radius':1.5})
        self.assertAlmostEqual(stats['mean'],.5);self.assertAlmostEqual(stats['median'],.5)
        self.assertEqual(stats['valid_samples'],8);self.assertEqual(stats['field_samples'],9)
    def test_dnb_conversion_uses_extended_band(self):
        w=np.arange(300.,901.);src=np.c_[w,np.ones(len(w))];rsr=np.array([[300,1.],[900,1.]])
        self.assertAlmostEqual(overlap_conversion(src,rsr,0,band=(352,750))[0],398/600)

if __name__=='__main__':unittest.main()
