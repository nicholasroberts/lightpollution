"""Regression for a wrong polar-rotation start, as observed in Gavarnie image 12."""
import sys
from pathlib import Path
import unittest
import numpy as np
from scipy import ndimage
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from analyse_spider_sky_rotation import refine_multistart

class MultistartTests(unittest.TestCase):
    def test_nominal_start_recovers_from_wrong_polar_peak(self):
        rng=np.random.default_rng(82)
        image=10+ndimage.gaussian_filter(rng.normal(size=(200,200)),1)
        moving=ndimage.shift(ndimage.rotate(image,-135,reshape=False,order=1),[5,-4],order=1)
        fit,_=refine_multistart(image,moving,[99.5]*2,[99.5]*2,65,174,135)
        self.assertLess(abs(fit['scene_roll_clockwise_deg']-135),.3)
        self.assertGreater(fit['inner_sky_correlation'],.9)
        self.assertLess(np.linalg.norm(np.array(fit['source_center_correction_xy_px'])-[-4,5]),1)
        self.assertEqual(len(fit['multistart_trials']),3)

if __name__=='__main__':unittest.main()
