import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wang2016_sky import (complex_pattern,neutral_points,turbidity_multiplier,
    natural_stokes,add_unpolarized,equivalent_turbidity_increment,spectral_factor,horizon_factor)

class WangTests(unittest.TestCase):
    def test_neutral_points_zero(self):
        a,b=neutral_points(.7,np.pi/2,24)
        for root in [a,b]:
            w=complex_pattern(np.array(2*np.arctan(abs(root))),np.array(np.angle(root)),.7,np.pi/2,24)
            self.assertLess(abs(w),1e-12)
    def test_printed_turbidity_exponent(self):
        self.assertAlmostEqual(turbidity_multiplier(3),np.exp(-3/4+.12))
    def test_horizon_endpoints(self):
        np.testing.assert_allclose(horizon_factor(np.array([0,np.pi/2])),[1,0],atol=.025)
    def test_alan_preserves_qu_and_dilutes_dop(self):
        i=np.array([1.,2.]);q=np.array([.3,.4]);u=np.array([.2,-.1])
        ii,qq,uu=add_unpolarized(i,q,u,2*i)
        np.testing.assert_array_equal(qq,q);np.testing.assert_array_equal(uu,u)
        np.testing.assert_allclose(np.hypot(qq,uu)/ii,np.hypot(q,u)/i/3)
        np.testing.assert_allclose(np.arctan2(uu,qq),np.arctan2(u,q))
    def test_equivalence_is_proportional_intensity_only(self):
        delta=equivalent_turbidity_increment(2)
        self.assertAlmostEqual(turbidity_multiplier(3+delta)/turbidity_multiplier(3),1/3)
    def test_spectral_factor_anchors(self):
        self.assertAlmostEqual(float(spectral_factor(2.,2.,4.)),1.)
        self.assertAlmostEqual(float(spectral_factor(4.,2.,4.)),0.)
    def test_full_equation_requires_spectral_input(self):
        theta=np.linspace(0,1,20);phi=np.zeros_like(theta)
        with self.assertRaises(ValueError):natural_stokes(theta,phi,.7,.4,24,1,mode='eq17')
        i,q,u=natural_stokes(theta,phi,.7,.4,24,1,mode='eq17',spectral_s=np.full_like(theta,.5))
        self.assertTrue(np.all(np.hypot(q,u)<=i))
    def test_invalid_alan_rejected(self):
        with self.assertRaises(ValueError):add_unpolarized(1,.2,.1,-1)

if __name__=='__main__':unittest.main()
