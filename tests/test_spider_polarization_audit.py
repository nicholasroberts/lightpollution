import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from audit_spider_polarization import summarize

class AuditTests(unittest.TestCase):
    def test_spatial_cancellation(self):
        i=np.ones((10,10));q=np.full_like(i,.4);q[:5]*=-1;u=np.zeros_like(i)
        s,_,_=summarize(i,q,u,np.ones_like(i,dtype=bool))
        self.assertAlmostEqual(s['physical_mean_dop'],.4)
        self.assertAlmostEqual(s['fixed_eye_results']['stokes_intensity']['prc'],0)
        self.assertGreater(s['hypothetical_locally_aligned_benchmark']['stokes_intensity_prc'],0)

    def test_axial_wrap_does_not_cancel(self):
        angle=np.deg2rad(np.array([[-89.,89.]]));i=np.ones_like(angle)
        s,_,_=summarize(i,.4*np.cos(2*angle),.4*np.sin(2*angle),np.ones_like(i,dtype=bool))
        self.assertGreater(s['orientation_coherence_unweighted'],.99)

    def test_invalid_values_reported(self):
        i=np.ones((2,3));q=np.array([[.2,2.,np.nan],[.2,.2,.2]]);u=np.zeros_like(i)
        s,_,_=summarize(i,q,u,np.ones_like(i,dtype=bool))
        self.assertEqual(s['dop_above_one_pixels'],1)
        self.assertEqual(s['nonfinite_or_nonpositive_pixels'],1)
        self.assertEqual(s['physical_pixels'],4)

if __name__=='__main__':unittest.main()
