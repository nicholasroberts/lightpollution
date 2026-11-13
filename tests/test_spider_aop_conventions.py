import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from audit_spider_aop_conventions import transform,fit_offset,residual_stats,apparent_eta

class ConventionTests(unittest.TestCase):
    def test_conjugation_and_common_offset(self):
        phi=np.linspace(-2,2,100);model=.6*np.exp(2j*phi)
        obs=np.conj(model)*np.exp(2j*.3)
        offset=fit_offset(obs,transform(model,phi,'conjugate'))
        self.assertAlmostEqual(offset,.3)
        other=.4*np.exp(2j*np.linspace(-1,2,55))
        heldout=np.conj(other)*np.exp(2j*.3)
        self.assertLess(residual_stats(heldout,np.conj(other)*np.exp(2j*offset))['mae_deg'],1e-10)
        self.assertGreater(residual_stats(heldout,other*np.exp(2j*offset))['mae_deg'],20)
    def test_basis_spin(self):
        phi=np.linspace(-2,2,15);z=np.ones(15,dtype=complex)
        np.testing.assert_allclose(transform(z,phi,'direct_radial_plus'),np.exp(2j*phi))
        np.testing.assert_allclose(abs(transform(z,phi,'conjugate_radial_minus')),1)
    def test_atmosphere_alan_separate(self):
        m=np.array([.3,.6,.9]);obs=.8*m/2
        unknown=apparent_eta(obs,m)
        self.assertNotIn('conditional_eta_quantiles',unknown)
        np.testing.assert_allclose(unknown['observed_over_natural_quantiles'],[.4]*3)
        np.testing.assert_allclose(apparent_eta(obs,m,1)['conditional_eta_quantiles'],[.8]*3)
    def test_reject_negative_alan(self):
        with self.assertRaises(ValueError):apparent_eta(np.array([.2]),np.array([.5]),-1)

if __name__=='__main__':unittest.main()
