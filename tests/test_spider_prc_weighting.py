import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from audit_spider_prc_weighting import analyze,weighted_pattern


class WeightingTests(unittest.TestCase):
    def test_factorization_and_cancellation(self):
        r=weighted_pattern(np.array([1.,1.]),np.array([.8,-.8],complex),np.ones(2,bool))
        self.assertAlmostEqual(r['weighted_mean_local_dolp'],.8)
        self.assertAlmostEqual(r['aop_cancellation_factor'],0)
        self.assertEqual(r['maximum_absolute_prc'],0)
        self.assertEqual(r['integrated_dolp_factorization_error'],0)

    def test_common_weights_remove_intensity_only_difference(self):
        z=np.array([.2+.1j,-.1+.2j]);i=np.array([1.,9.]);j=i[::-1];mask=np.ones(2,bool)
        matrix,_=analyze(dict(measured=(i,(i*z).real,(i*z).imag),mixed=(j,(j*z).real,(j*z).imag)),mask)
        for weight in matrix:
            self.assertAlmostEqual(matrix[weight]['measured']['maximum_absolute_prc'],matrix[weight]['mixed']['maximum_absolute_prc'])
        self.assertNotAlmostEqual(matrix['measured']['measured']['sum_Q'],matrix['mixed']['mixed']['sum_Q'])

    def test_zero_and_weight_scale_invariance(self):
        w=np.array([1.,2.,3.]);z=np.array([.2+.1j,.3-.2j,.1j]);mask=np.ones(3,bool)
        r=weighted_pattern(w,z,mask);s=weighted_pattern(7*w,z*np.exp(1.2j),mask)
        for key in ['integrated_dolp','weighted_mean_local_dolp','aop_cancellation_factor','maximum_absolute_prc']:
            self.assertAlmostEqual(r[key],s[key])
        self.assertLess(abs(r['integrated_dolp_factorization_error']),1e-14)


if __name__=='__main__':unittest.main()
