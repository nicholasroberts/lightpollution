"""Checks coordinate signs, rotation recovery and pipeline export on synthetic skies."""
import sys
from pathlib import Path
import tempfile
import json
import unittest
from unittest.mock import patch
import numpy as np
import tifffile
from scipy import ndimage
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import analyse_spider_sky_rotation as rotation
from reconstruct_spider_eye_contrast import digest

class RotationTests(unittest.TestCase):
    def test_fft_sign_and_global_peaks(self):
        rng=np.random.default_rng(23); ref=rng.normal(size=(200,720))
        for shift in [0,90,-90,180,270,-270]:
            fit=rotation.estimate_roll(ref,np.roll(ref,shift,axis=1))
            error=(fit['scene_roll_clockwise_deg']-shift/2+180)%360-180
            self.assertLess(abs(error),.01)

    def test_coordinate_warp(self):
        y,x=np.indices((101,101)); im=np.exp(-((x-70)**2+(y-50)**2)/10)
        aligned=rotation.warp_to_reference(im,[50,50],[50,50],90)
        self.assertEqual(np.unravel_index(np.nanargmax(aligned),im.shape),(30,50))

    def test_scene_structure_recovers_roll(self):
        rng=np.random.default_rng(9); im=ndimage.gaussian_filter(rng.normal(size=(256,256)),1)+10
        im*=np.hypot(*np.indices(im.shape)-127.5)<115
        ref=rotation.polar_features(im,[127.5,127.5],115)
        for angle in [45,90,135,-45]:
            moved=ndimage.rotate(im,-angle,reshape=False,order=1)
            found=rotation.estimate_roll(ref,rotation.polar_features(moved,[127.5,127.5],115))
            self.assertLess(abs(found['scene_roll_clockwise_deg']-angle),.2)

    def test_single_quartet_end_to_end(self):
        rng=np.random.default_rng(42)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); folder=root/'inputs';folder.mkdir()
            y,x=np.indices((160,160)); sky=100+ndimage.gaussian_filter(rng.normal(size=x.shape)*20,1)
            sky[np.hypot(x-79.5,y-79.5)>73]=0
            records=[]
            for num,angle in enumerate([0,45,90,135],1):
                # Uniform polarization with known modulation; geometry structural match independent of gain.
                im=ndimage.rotate(sky*(1+.2*np.cos(np.deg2rad(2*angle))),-angle,reshape=False,order=1)
                path=folder/f'{num}.tif';tifffile.imwrite(path,np.stack([im]*3,axis=-1).astype('float32'))
                records.append(dict(file=str(path),sha256=digest(path),fits=[dict(threshold_fraction=.05,
                    center_x_px=79.5,center_y_px=79.5,rim_radius_px=73)]))
            geometry=root/'geometry.json';geometry.write_text(json.dumps(dict(channel=2,records=records)))
            args=['rotation','--folder',str(folder),'--geometry',str(geometry),'--output',str(root/'out')]
            with patch.object(sys,'argv',args),patch('reconstruct_spider_eye_contrast.reference_comparison',return_value={'status':'synthetic'}):
                rotation.main()
            report=json.loads((root/'out/rotation_analysis.json').read_text())
            self.assertEqual(len(report['rotations']),4)
            self.assertEqual(len(report['prc_diagnostics']),12)
            with np.load(root/'out/linear_stokes.npz') as stokes:
                mask=(x-79.5)**2+(y-79.5)**2<45**2
                self.assertAlmostEqual(np.mean(stokes['Q'][mask]/stokes['I'][mask]),.2,delta=.01)
                self.assertLess(abs(np.mean(stokes['U'][mask]/stokes['I'][mask])),.01)
            self.assertTrue((root/'out/linear_stokes_preview/result.json').exists())
            self.assertTrue((root/'out/legacy_dop_aop_preview/result.json').exists())

if __name__=='__main__': unittest.main()
