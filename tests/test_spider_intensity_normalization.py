import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from audit_spider_intensity_normalization import ideal_rayleigh,mixed_dop
from wang2016_sky import berry_normalized_pattern

class NormalizationTests(unittest.TestCase):
    def test_rayleigh_half_maximum(self):
        angle=np.arccos(1/np.sqrt(3))
        self.assertAlmostEqual(float(ideal_rayleigh(angle)),.5)
        self.assertAlmostEqual(float(ideal_rayleigh(np.pi-angle)),.5)
    def test_berry_coincident_limit(self):
        g=np.linspace(.01,np.pi-.01,501)
        w=berry_normalized_pattern(np.full_like(g,np.pi/2),np.pi/2+g,np.pi/2,np.pi/2,.001)
        np.testing.assert_allclose(abs(w),np.sin(g)**2,atol=1e-9)
    def test_intensity_rescale_invariance(self):
        i=np.array([1,2,3]);p=np.array([.2,.5,.7]);a=.4
        np.testing.assert_allclose(mixed_dop(i,p,a),mixed_dop(i*5,p,a*5))
    def test_unpolarized_addition(self):
        self.assertAlmostEqual(float(mixed_dop(2,.6,2)),.3)
        with self.assertRaises(ValueError):mixed_dop(1,1.2,0)

if __name__=='__main__':unittest.main()
