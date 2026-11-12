import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wang2016_sky import natural_stokes,add_unpolarized,alan_relative_from_depth

class AlanOnlyTests(unittest.TestCase):
    def test_default_has_no_atmospheric_multiplier(self):
        args=(np.array([.2,.4]),np.array([.1,.3]),1.,1.,24.,np.ones(2))
        default=natural_stokes(*args)
        explicit=natural_stokes(*args,depolarization_multiplier=1.)
        for a,b in zip(default,explicit):np.testing.assert_allclose(a,b)
    def test_depth_is_exact_additive_dilution(self):
        i=np.array([1.,3.]);q=np.array([.2,.4]);u=np.array([.1,-.2])
        depth=.7;alpha=alan_relative_from_depth(depth)
        ii,qq,uu=add_unpolarized(i,q,u,alpha*i)
        np.testing.assert_allclose((qq+1j*uu)/ii,(q+1j*u)/i*np.exp(-depth))
        np.testing.assert_array_equal(qq,q);np.testing.assert_array_equal(uu,u)
        self.assertEqual(float(alan_relative_from_depth(0)),0.)
    def test_negative_depth_rejected(self):
        with self.assertRaises(ValueError):alan_relative_from_depth(-1.)
