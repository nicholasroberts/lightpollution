import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from fit_spider_prc_lamp_equivalence import required_energy


class LampEquivalenceTests(unittest.TestCase):
    def test_known_prc_target(self):
        b,g=.3,.2;ni,p,c,n=20.,4.,2.,10
        amount=3.;target=2*np.arctanh(g*p/(b*(ni+amount*c*n)))
        fit=required_energy(target,ni,p,b,g,c,n)
        self.assertAlmostEqual(fit['energy'],amount)
        self.assertAlmostEqual(fit['achieved_prc'],target)

    def test_higher_receptor_overlap_requires_less_energy(self):
        h=required_energy(.1,20,4,.3,.2,1,10)
        l=required_energy(.1,20,4,.3,.2,2,10)
        self.assertAlmostEqual(l['energy']/h['energy'],.5)

    def test_unreachable_target_and_zero_addition(self):
        natural=2*np.arctanh(.2*4/(.3*20))
        self.assertEqual(required_energy(natural+1,20,4,.3,.2,1,10)['status'],'unreachable_by_unpolarized_addition')
        self.assertAlmostEqual(required_energy(natural,20,4,.3,.2,1,10)['energy'],0)


if __name__=='__main__':unittest.main()
