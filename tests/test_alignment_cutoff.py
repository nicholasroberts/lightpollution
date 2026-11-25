"""A declared validation cutoff changes acceptance, not the geometric fit."""
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from spider_fisheye_registration import fit_correspondences


class AlignmentCutoffTests(unittest.TestCase):
    def test_explicit_four_pixels_changes_acceptance_only(self):
        angle=np.linspace(0,2*np.pi,64,endpoint=False)
        source=np.c_[100+80*np.cos(angle),100+80*np.sin(angle)]
        target=source.copy();target[::4,0]+=3.5
        center=np.array([100.,100.])
        old,q,e=fit_correspondences(source,target,center,100.)
        new,r,f=fit_correspondences(source,target,center,100.,
            validation_cutoff_px=4,validation_median_limit_px=4)
        self.assertFalse(q['numerical_pass'])
        self.assertTrue(r['numerical_pass'])
        np.testing.assert_allclose(e,f)
        np.testing.assert_allclose(old['ray_rotation_source_to_reference'],new['ray_rotation_source_to_reference'])
        np.testing.assert_allclose(old['target_projection_center_xy_px'],new['target_projection_center_xy_px'])
        self.assertEqual(r['validation_inlier_cutoff_px'],4)


if __name__=='__main__':unittest.main()
