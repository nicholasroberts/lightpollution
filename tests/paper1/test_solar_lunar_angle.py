"""Checks of the angle-structure mixture objective, independent of site data."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from fit_paper1_solar_lunar_angle import optimize_fraction


class MixtureTests(unittest.TestCase):
    def recover(self,fraction,offset):
        rng=np.random.default_rng(42)
        ml,mh,sl,sh=[rng.normal(size=400)+1j*rng.normal(size=400) for _ in range(4)]
        lo=(1-fraction)*ml+fraction*sl
        hi=(1-fraction)*mh+fraction*sh
        phase=np.exp(2j*offset)
        change=(hi*phase)*np.conj(lo*phase)
        fitted,*_=optimize_fraction(change,ml,mh,sl,sh,np.ones(400,dtype=bool))
        self.assertAlmostEqual(fitted,fraction,places=6)

    def test_mixed_with_unknown_zero(self):
        self.recover(.37,.7)

    def test_lunar_endpoint(self):
        self.recover(0.,-.8)

    def test_solar_endpoint(self):
        self.recover(1.,1.2)


if __name__=='__main__':
    unittest.main()
