import unittest,sys
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from spider_fisheye_registration import rays,pixels,forward,fit_correspondences,warp


class FisheyeRegistrationTest(unittest.TestCase):
    def test_projection_round_trip(self):
        points=np.array([[300.,300.],[450,300],[300,80],[200,470]])
        np.testing.assert_allclose(pixels(rays(points,np.array([300,300]),280),np.array([300,300]),280),points,atol=1e-10)

    def test_tilt_mapping_recovers_withheld_points_despite_outliers(self):
        rng=np.random.default_rng(12);angle=rng.uniform(-np.pi,np.pi,300);rad=rng.uniform(180,270,300)
        source=np.array([300.,300.])+np.column_stack([rad*np.cos(angle),rad*np.sin(angle)])
        mapping=dict(source_center_xy_px=[302.,298.],target_projection_center_xy_px=[299.,303.],radius_px=280.,
            ray_rotation_source_to_reference=Rotation.from_euler('xyz',[2,-1,45],degrees=True).as_matrix().tolist())
        target=forward(source,mapping)+rng.normal(0,.05,source.shape)
        target[1::9]=rng.uniform(50,550,(len(target[1::9]),2))
        fitted,q,error=fit_correspondences(source,target,np.array([300.,300.]),280.)
        self.assertTrue(q['numerical_pass'],q)
        self.assertLess(q['heldout_median_all_matches_px'],.4)
        self.assertAlmostEqual(np.linalg.det(fitted['ray_rotation_source_to_reference']),1.,places=8)

    def test_inverse_warp_roll_has_correct_direction(self):
        mapping=dict(source_center_xy_px=[50.,50.],target_projection_center_xy_px=[50.,50.],radius_px=45.,
            ray_rotation_source_to_reference=Rotation.from_euler('z',90,degrees=True).as_matrix().tolist())
        image=np.zeros((101,101));image[50,70]=1
        out=warp(image,mapping)
        self.assertAlmostEqual(out[70,50],1.,places=8)
        self.assertTrue(np.isnan(out[0,0]))


if __name__=='__main__':unittest.main()
