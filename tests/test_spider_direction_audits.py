import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from audit_spider_mixed_aop import common_mask,metrics
from audit_spider_analyzer_intensities import forward,audit


class DirectionTests(unittest.TestCase):
    def test_unknown_global_zero_cancels(self):
        y,x=np.indices((12,12));phase=np.exp(2j*(.02*x+.03*y))
        shifted=phase*np.exp(2j*1.2)
        result=metrics(shifted,phase,np.ones((12,12),bool),2)
        self.assertLess(result['p95_absolute_deg'],1e-10)

    def test_conjugation_changes_gradient(self):
        y,x=np.indices((12,12));phase=np.exp(2j*(.02*x+.03*y));mask=np.ones((12,12),bool)
        self.assertGreater(metrics(np.conj(phase),phase,mask,2)['median_absolute_deg'],1)
        self.assertLess(metrics(np.conj(np.conj(phase)),phase,mask,2)['median_absolute_deg'],1e-10)

    def test_identical_pairs_on_common_support(self):
        phase=np.ones((12,12),complex);phases={k:phase.copy() for k in ['moon','sun','mixed']}
        dolps={k:np.full((12,12),.2) for k in phases};dolps['moon'][2:6,2:6]=.01
        phases['sun'][7,7]=np.nan
        mask=common_mask(np.ones((12,12),bool),np.full((12,12),.2),phases,dolps)
        self.assertFalse(mask[3,3]);self.assertFalse(mask[7,7])
        counts=[metrics(phase,p,mask,2)['pair_count'] for p in phases.values()]
        self.assertEqual(len(set(counts)),1)

    def test_linear_analyzer_closure_and_known_stokes(self):
        stokes=np.stack([np.full((8,8),100.),np.full((8,8),20.),np.full((8,8),10.)])
        angles=[0,45,90,135];data=forward(stokes,angles)
        recovered,loo,summary=audit(data,angles,np.ones((8,8),bool))
        np.testing.assert_allclose(recovered,stokes,atol=1e-10)
        self.assertLess(summary['leave_one_out_fraction']['p95_absolute'],1e-12)


if __name__=='__main__':unittest.main()
