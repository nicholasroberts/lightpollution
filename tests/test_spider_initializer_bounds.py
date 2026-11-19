import sys
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import spider_shared_fisheye_registration as s

class InitializerTests(unittest.TestCase):
    def test_roundoff_outside_bounds(self):
        lo=np.array([-np.inf,-5.,-5.]);hi=np.array([np.inf,5.,5.])
        x=s.feasible_start([3.,5.+1e-12,-5.-1e-12],lo,hi)
        self.assertEqual(x[0],3.)
        self.assertTrue(np.all(x[1:]>lo[1:]) and np.all(x[1:]<hi[1:]))
        np.testing.assert_array_equal(s.feasible_start([3,1,2],lo,hi),[3,1,2])
        with self.assertRaises(ValueError):s.feasible_start([np.nan,1,2],lo,hi)

    def test_shared_fit_from_boundary_initializer(self):
        rng=np.random.default_rng(123);center=np.array([1000.,666.]);radius=650.
        theta=rng.uniform(-np.pi,np.pi,250);rr=rng.uniform(.6,.98,250)*radius
        points=center+np.c_[rr*np.cos(theta),rr*np.sin(theta)]
        mapping={'ray_rotation_source_to_reference':np.eye(3).tolist(),
                 'target_projection_center_xy_px':(center+[.05*radius+1e-12,0]).tolist()}
        with patch.object(s,'fit_correspondences',return_value=(mapping,{},None)):
            fits,shared=s.fit_site({4:(points,points),7:(points,points)},center,radius)
        self.assertTrue(all(x[1]['numerical_pass'] for x in fits.values()))
        self.assertGreater(shared['joint_initializer_adjustment_max'],0)

if __name__=='__main__':unittest.main()
