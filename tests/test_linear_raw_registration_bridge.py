import json,sys,tempfile,unittest,hashlib
from pathlib import Path
from unittest.mock import patch
import numpy as np
from PIL import Image
import tifffile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lock_spider_site_registration import process_site


class LinearRawBridgeTest(unittest.TestCase):
    def test_missing_seven_uses_two_complete_quartets_and_new_canvas(self):
        with tempfile.TemporaryDirectory() as root:
            root=Path(root);batch=root/'batch';old_dir=root/'old';raw_dir=root/'linear';out=root/'out'
            for p in [old_dir,raw_dir,out,batch/'museum/registration']:p.mkdir(parents=True)
            inputs=[];frames=[];shutter='2019:05:14 22:09:51'
            values={2:10,5:20,8:30,11:40,3:50,6:60,9:70,12:80}
            def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
            for n in range(1,13):
                p=old_dir/f'{n}.tif'
                Image.fromarray(np.zeros((133,200,3),np.uint8)).save(p)
                inputs.append(dict(path=str(p),sha256=digest(p),shape=[133,200,3],channel=2))
                if n==7:continue
                p=raw_dir/f'{n}.tif';tifffile.imwrite(p,np.full((132,200),values.get(n,1000),np.float32),metadata=None)
                frames.append(dict(exposure=n,datetime_original=shutter,source_sha256='synthetic_raw_fixture',output_sha256=digest(p)))
            (batch/'museum/registration/rotation_analysis.json').write_text(json.dumps(dict(
                site_id='test',reference_center_xy_px=[100,66],radius_60_px=65/np.sqrt(2),inputs=inputs)))
            manifest=raw_dir/'linear_raw_manifest.json';manifest.write_text(json.dumps(dict(
                status='linear_raw_blue_prepared_not_registered_complete',frames=frames,
                processing=dict(method='visible_Bayer_blue_samples_minus_camera_reported_black',
                    white_balance_applied=False,gamma_applied=False,color_matrix_applied=False,
                    display_profile_or_look_applied=False,auto_brightness_applied=False),
                excluded_exposures=[7],analyzer_groups=[[2,5,8,11],[3,6,9,12]])))
            points=np.array([[40,30],[100,20],[160,30],[170,65],[160,100],[100,110],[40,100],[30,65]])
            def fit(pairs,c,r):
                fitted={n:(dict(source_center_xy_px=c.tolist(),target_projection_center_xy_px=c.tolist(),radius_px=r,
                    ray_rotation_source_to_reference=np.eye(3).tolist(),shared_tangential_p1_p2=[0,0]),
                    dict(numerical_pass=True),np.zeros(len(s))) for n,(s,t) in pairs.items()}
                self.assertNotIn(7,pairs)
                return fitted,{'fixture_only':True}
            with patch('lock_spider_site_registration.landmark_features',return_value=None),\
                 patch('lock_spider_site_registration.match_landmarks',return_value=(points,points)),\
                 patch('spider_shared_fisheye_registration.fit_site',side_effect=fit),\
                 patch('lock_spider_site_registration.moon_coordinate',return_value=(None,{})),\
                 patch('PIL.Image.open') as opened:
                opened.return_value.__enter__.return_value.getexif.return_value.get_ifd.return_value={36867:shutter}
                lock=process_site('museum',batch,out,None,{},registration_model='fisheye-shared',linear_raw_manifest=manifest)
            self.assertEqual(lock['analyzer_groups'],[[2,5,8,11],[3,6,9,12]])
            self.assertEqual(lock['excluded_exposures'],[7]);self.assertEqual(len(lock['transforms']),11)
            self.assertEqual(lock['geometry_transfer']['new_shape'],[132,200])
            self.assertEqual(lock['radiometry_status'],'black_subtracted_linear_camera_blue')
            arr=np.load(out/'museum/registration/linear_stokes.npz')
            self.assertAlmostEqual(arr['I'][66,100],80.,places=5)
            self.assertAlmostEqual(arr['Q'][66,100],-20.,places=5)
            self.assertAlmostEqual(arr['U'][66,100],-20.,places=5)


if __name__=='__main__':unittest.main()
