import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from diagnose_spider_replicates import replicate_results

class ReplicateDiagnosticTests(unittest.TestCase):
    def quartet(self,sign=1):
        return [np.full((8,8),v) for v in [60 if sign==1 else 40,50,40 if sign==1 else 60,50]]
    def test_identical_replicates_preserve_contrast(self):
        data=replicate_results([self.quartet()]*3,np.ones((8,8),dtype=bool))
        for entries in data.values():
            self.assertEqual(len(entries),4)
            values=[e['fixed_eye_results']['stokes_intensity']['prc'] for e in entries]
            self.assertTrue(np.allclose(values,values[0]))
    def test_distinct_orientations_cancel_only_when_combined(self):
        data=replicate_results([self.quartet(),self.quartet(-1)],np.ones((8,8),dtype=bool))
        entries=data['linear_stokes']
        self.assertGreater(entries[0]['fixed_eye_results']['stokes_intensity']['prc'],.2)
        self.assertAlmostEqual(entries[-1]['fixed_eye_results']['stokes_intensity']['prc'],0)

if __name__=='__main__':unittest.main()
