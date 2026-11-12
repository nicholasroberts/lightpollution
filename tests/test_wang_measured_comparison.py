import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from compare_wang2016_measured_skies import axial_offset_and_residual,source_candidate

class MeasuredComparisonTests(unittest.TestCase):
    def test_axial_offset_preserves_amplitude(self):
        angle=np.linspace(-1,1,100).reshape(10,10)
        model=.4*np.exp(2j*angle);observed=.2*np.exp(2j*(angle+.3))
        offset,res=axial_offset_and_residual(observed,model,np.ones(angle.shape,dtype=bool))
        self.assertAlmostEqual(offset,.3)
        self.assertLess(res['median_abs_deg'],1e-10)
        np.testing.assert_allclose(abs(model),.4)
    def test_compact_source_detection(self):
        y,x=np.indices((100,100));image=10+100*np.exp(-((x-65)**2+(y-40)**2)/8)
        position,_=source_candidate(image,np.ones(image.shape,dtype=bool))
        np.testing.assert_allclose(position,[65,40])

if __name__=='__main__':unittest.main()
