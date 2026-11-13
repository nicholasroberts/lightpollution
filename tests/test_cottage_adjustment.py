import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from adjust_cottage_alan import target_addition,mean_dop

class CottageAdjustmentTests(unittest.TestCase):
    def test_known_target(self):
        pn=np.full(100,.675);i=np.ones(100)
        alpha=target_addition(pn,i,.45)
        self.assertAlmostEqual(alpha,.5)
        self.assertAlmostEqual(mean_dop(pn,i,alpha,1,1),.45)
    def test_stronger_levels_reduce_dop(self):
        pn=np.linspace(.5,.8,100);i=np.linspace(.8,1.2,100)
        alpha=target_addition(pn,i,.45,2.,3.)
        values=[mean_dop(pn,i,alpha*level/5,2.,3.) for level in [5,6,7]]
        self.assertAlmostEqual(values[0],.45)
        self.assertGreater(values[0],values[1]);self.assertGreater(values[1],values[2])
    def test_unreachable_target_rejected(self):
        with self.assertRaises(ValueError):target_addition(np.full(10,.3),np.ones(10),.45)

if __name__=='__main__':unittest.main()
