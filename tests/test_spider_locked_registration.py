import sys,unittest,tempfile,json,importlib.util
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lock_spider_site_registration import rigid_fit,robust_rigid,transform_points,warp,north_rotation,landmark_features,register_landmarks,moon_coordinate
from plot_spider_site_details import locked_measured
from reconstruct_spider_eye_contrast import digest


class LockedRegistrationTest(unittest.TestCase):
    def test_rigid_recovery_with_false_matches(self):
        rng=np.random.default_rng(2);source=rng.uniform(0,500,(100,2))
        angle=np.deg2rad(47);rotation=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        truth=np.column_stack([rotation,[13.,-7.]])
        target=transform_points(source,truth)+rng.normal(0,.1,(100,2))
        target[:30]=rng.uniform(0,500,(30,2))
        fit,mask=robust_rigid(source,target)
        np.testing.assert_allclose(fit,truth,atol=.08)
        self.assertGreaterEqual(mask.sum(),68)
        self.assertAlmostEqual(np.linalg.det(fit[:,:2]),1)

    def test_warp_and_north_do_not_mirror(self):
        center=np.array([50.,50.]);moon=np.array([80.,50.])
        matrix,angle=north_rotation(center,moon,0.)
        np.testing.assert_allclose(transform_points([moon],matrix),[[50,20]],atol=1e-10)
        self.assertAlmostEqual(np.linalg.det(matrix[:,:2]),1.)
        im=np.zeros((101,101));im[50,80]=1;im[70,50]=.5
        out=warp(im,matrix)
        self.assertAlmostEqual(out[20,50],1.)
        self.assertAlmostEqual(out[50,70],.5)

    def test_moon_selection_independent_of_ephemeris(self):
        y,x=np.indices((201,201));im=np.ones((201,201))+100*np.exp(-((x-170)**2+(y-95)**2)/8)
        center=np.array([100.,100.])
        a,ma=moon_coordinate(im,center,100.,'minifactory',{'moon_elevation_deg':5},{})
        b,mb=moon_coordinate(im,center,100.,'minifactory',{'moon_elevation_deg':70},{})
        np.testing.assert_equal(a,b)
        self.assertNotEqual(ma['elevation_residual_deg'],mb['elevation_residual_deg'])

    @unittest.skipUnless(importlib.util.find_spec('cv2'), 'Optional registration dependency not installed')
    def test_real_feature_pipeline_on_rotated_static_texture(self):
        rng=np.random.default_rng(8);shape=(601,601);center=np.array([300.,300.]);radius=280.
        texture=gaussian_filter(rng.normal(0,1,shape),1)
        ref=100+30*texture
        for _ in range(200):
            px,py=rng.integers(30,570,2);ref[py-2:py+3,px-2:px+3]+=rng.uniform(10,100)
        a=np.deg2rad(44);rotation=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
        ref_to_moving=np.column_stack([rotation,center-rotation@center+[3,-4]])
        moving=warp(ref,ref_to_moving)
        matrix,q,_,_,_=register_landmarks(landmark_features(ref,center,radius),landmark_features(moving,center,radius),center,radius)
        self.assertTrue(q['numerical_pass'],q)
        expected=np.column_stack([rotation.T,-rotation.T@ref_to_moving[:,2]])
        np.testing.assert_allclose(matrix,expected,atol=.3)

    def test_frozen_reader_rejects_changed_products_and_astronomy(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);folder=root/'site/registration';folder.mkdir(parents=True)
            product=folder/'north_up_dop.npz';np.savez(product,DoLP=np.full((101,101),.3),I=np.ones((101,101)))
            anchor=dict(utc='2019-01-01T22:00:00',timezone='Europe/London',moon_elevation_deg=20.,moon_compass_azimuth_deg=180.)
            record=dict(status='numerically_locked_awaiting_visual_review',astronomical_anchor=anchor,
                reference_center_xy_px=[50,50],rim_radius_px=45,north_up_dop_path='site/registration/north_up_dop.npz')
            lock=root/'locked_site_references.json';lock.write_text(json.dumps(dict(sites=dict(site=record))))
            checks={lock.name:digest(lock),str(product.relative_to(root)):digest(product)}
            (root/'reference_sha256.json').write_text(json.dumps(checks))
            arr,meta,north=locked_measured(lock,'site',anchor,grid=41)
            self.assertTrue(north);np.testing.assert_allclose(arr[np.isfinite(arr)],.3)
            with self.assertRaisesRegex(ValueError,'anchor differs'):
                locked_measured(lock,'site',dict(anchor,moon_compass_azimuth_deg=90))
            np.savez(product,DoLP=np.full((101,101),.8),I=np.ones((101,101)))
            with self.assertRaisesRegex(ValueError,'hash mismatch'):locked_measured(lock,'site',anchor)


if __name__=='__main__':unittest.main()
