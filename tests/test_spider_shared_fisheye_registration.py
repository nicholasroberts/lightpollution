import sys, unittest
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from spider_shared_fisheye_registration import undistort,distort,forward,inverse_points,warp,heldout_cells,fit_site


class SharedFisheyeTest(unittest.TestCase):
    def test_inverse_and_warp_direction(self):
        center=np.array([50.,50.]);radius=45.;k=[.002,.001]
        points=np.array([[50.,50.],[70,50],[25,65],[42,12]])
        np.testing.assert_allclose(distort(undistort(points,center,radius,k),center,radius,k),points,atol=1e-9)
        mapping=dict(source_center_xy_px=center.tolist(),target_projection_center_xy_px=center.tolist(),radius_px=radius,
            shared_tangential_p1_p2=k,ray_rotation_source_to_reference=Rotation.from_euler('z',90,degrees=True).as_matrix().tolist())
        np.testing.assert_allclose(inverse_points(forward(points,mapping),mapping),points,atol=1e-9)
        identity=dict(mapping,ray_rotation_source_to_reference=np.eye(3).tolist())
        image=np.arange(101*101).reshape(101,101).astype(float);result=warp(image,identity)
        np.testing.assert_allclose(result[30:70,30:70],image[30:70,30:70],atol=1e-8)
        self.assertTrue(np.isnan(result[0,0]))

    def test_shared_fit_transfers_to_untrained_rotation(self):
        rng=np.random.default_rng(42);center=np.array([1000.,660.]);radius=650.;pairs={};k=[.002,.001]
        for n,angle in [(4,45),(7,90),(10,135)]:
            theta=rng.uniform(-np.pi,np.pi,360);rad=rng.uniform(.55,.97,360)*radius
            source=center+np.c_[rad*np.cos(theta),rad*np.sin(theta)]
            mapping=dict(source_center_xy_px=center.tolist(),target_projection_center_xy_px=(center+[3,4]).tolist(),radius_px=radius,
                shared_tangential_p1_p2=k,ray_rotation_source_to_reference=Rotation.from_euler('xyz',[2,-1,angle],degrees=True).as_matrix().tolist())
            target=forward(source,mapping)+rng.normal(0,.04,source.shape)
            pairs[n]=(source,target)
        fitted,shared=fit_site(pairs,center,radius)
        self.assertEqual(shared['calibration_exposures'],[4,7])
        self.assertTrue(shared['numerical_geometry_check'],shared)
        np.testing.assert_allclose(shared['coefficients_p1_p2'],k,atol=.0002)
        for mapping,q,error in fitted.values():
            self.assertTrue(q['numerical_pass'],q)
            self.assertLess(q['heldout_median_all_matches_px'],.2)
            self.assertGreater(np.linalg.det(mapping['ray_rotation_source_to_reference']),.999999)
        self.assertFalse(fitted[10][1]['shared_geometry_training_exposure'])
        # Entire cells have a single split membership, even for near-duplicates.
        self.assertEqual(heldout_cells([[101,101],[102,102]])[0],heldout_cells([[101,101],[102,102]])[1])


if __name__=='__main__':unittest.main()
