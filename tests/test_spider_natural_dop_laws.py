import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from compare_spider_natural_dop_laws import rayleigh_normalized_candidate,matched_profiles

class LawComparisonTests(unittest.TestCase):
    def test_rayleigh_limit(self):
        g=np.linspace(0,np.pi,100)
        np.testing.assert_allclose(abs(rayleigh_normalized_candidate(np.sin(g)**2,g)),np.sin(g)**2/(1+np.cos(g)**2))
    def test_phase_and_zeros_unchanged(self):
        w=np.array([0,.3+.2j,-.5j]);out=rayleigh_normalized_candidate(w,[0,.5,1.])
        self.assertEqual(out[0],0)
        np.testing.assert_allclose(np.angle(out[1:]),np.angle(w[1:]))
    def test_all_laws_share_pixels(self):
        g=np.full(30,45.);valid=np.ones(30,dtype=bool);valid[:5]=False
        obs=np.arange(30)/30;current=np.full(30,.8);candidate=np.full(30,.4)
        r=matched_profiles(g,obs,current,candidate,valid)[0]
        self.assertEqual(r['pixels'],25)
        self.assertAlmostEqual(r['measured_dop'],np.mean(obs[5:]))
        self.assertAlmostEqual(r['candidate_dop'],.4)

if __name__=='__main__':unittest.main()
