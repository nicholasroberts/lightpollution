import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from compare_spider_mixed_prc import integrated_prc,curve


class IntegratedPRCTest(unittest.TestCase):
    def test_cancellation_before_contrast(self):
        i=np.ones((1,2));q=np.array([[.8,-.8]]);u=np.zeros_like(i)
        r=integrated_prc(i,q,u,np.ones_like(i,dtype=bool))
        self.assertEqual(r['maximum_absolute_prc'],0)

    def test_intensity_weighted_not_mean_dolp(self):
        i=np.array([[1.,9.]]);q=np.array([[1.,0.]]);u=np.zeros_like(i)
        r=integrated_prc(i,q,u,np.ones_like(i,dtype=bool))
        self.assertAlmostEqual(r['integrated_dolp'],.1)
        self.assertNotAlmostEqual(r['integrated_dolp'],np.mean(q/i))

    def test_zero_rotation_and_scale_invariance(self):
        i=np.array([[4.,5.]]);z=np.array([[1+2j,2-1j]]);mask=np.ones_like(i,dtype=bool)
        r=integrated_prc(i,z.real,z.imag,mask)
        phase=np.exp(2j*np.deg2rad(37));rot=z*phase
        s=integrated_prc(8*i,8*rot.real,8*rot.imag,mask)
        for key in ['maximum_absolute_prc','continuous_heading_median_absolute_prc','integrated_dolp']:
            self.assertAlmostEqual(r[key],s[key])
        self.assertAlmostEqual((s['maximum_axis_native_deg']-r['maximum_axis_native_deg'])%180,37)
        dense=abs(curve(r,np.linspace(0,180,1000000,endpoint=False)))
        self.assertAlmostEqual(np.median(dense),r['continuous_heading_median_absolute_prc'],places=5)

    def test_matches_direct_receptor_absorption(self):
        i=np.array([[2.,4.]]);q=np.array([[.3,1.]]);u=np.array([[.5,-.2]])
        r=integrated_prc(i,q,u,np.ones_like(i,dtype=bool))
        theta=np.deg2rad(23)
        b=1-(10**(-.4)+10**(-.04))/2;a=(10**(-.04)-10**(-.4))/2
        first=np.sum(b*i+a*(q*np.cos(2*theta)+u*np.sin(2*theta)))
        second=np.sum(b*i-a*(q*np.cos(2*theta)+u*np.sin(2*theta)))
        self.assertAlmostEqual(float(curve(r,23)),np.log(first/second))

    def test_nonphysical_support_rejected(self):
        with self.assertRaises(ValueError):integrated_prc(np.ones((1,1)),np.ones((1,1))*2,np.zeros((1,1)),np.ones((1,1),bool))


if __name__=='__main__':unittest.main()
