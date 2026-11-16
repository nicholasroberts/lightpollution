import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from audit_spider_analyzer_intensities import solve,forward,leave_one_out,audit

class AnalyzerTest(unittest.TestCase):
    def test_nonorthogonal_angles_and_omitted_prediction(self):
        y,x=np.indices((12,14));s=np.array([10+0*x,2+.02*x,1+.03*y])
        angles=[0,42.77,86.07,132.77];d=forward(s,angles)
        np.testing.assert_allclose(solve(d,angles),s,atol=1e-12)
        np.testing.assert_allclose(leave_one_out(d,angles),d,atol=1e-12)
        mask=np.ones(x.shape,dtype=bool)
        _,_,good=audit(d,angles,mask);_,_,nominal=audit(d,[0,45,90,135],mask)
        self.assertLess(good['leave_one_out_fraction']['p95_absolute'],1e-12)
        self.assertGreater(nominal['leave_one_out_fraction']['median_absolute'],.001)
    def test_corrupted_frame_is_detected_and_sign_is_degenerate(self):
        s=np.array([np.ones((4,5))*10,np.ones((4,5))*2,np.ones((4,5))])
        angles=np.array([0,45,90,135]);d=forward(s,angles)
        negative=solve(d,-angles)
        np.testing.assert_allclose(negative[0:2],s[0:2],atol=1e-12)
        np.testing.assert_allclose(negative[2],-s[2],atol=1e-12)
        d[2]*=1.1
        self.assertGreater(np.max(abs(d-leave_one_out(d,angles))),.1)

if __name__=='__main__':unittest.main()
