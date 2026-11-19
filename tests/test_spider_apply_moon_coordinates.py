import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import apply_spider_moon_coordinates as m


class MoonCoordinateTests(unittest.TestCase):
    def test_geometry_preserves_center_and_fixed_bearing(self):
        matrix,meta=m.geometry([50,50],40,[70,50],{'moon_elevation_deg':40,'moon_compass_azimuth_deg':210})
        np.testing.assert_allclose(matrix[:,:2]@[50,50]+matrix[:,2],[50,50])
        self.assertAlmostEqual(np.linalg.det(matrix[:,:2]),1)
        v=matrix[:,:2]@[20,0]
        self.assertAlmostEqual((np.degrees(np.arctan2(v[1],v[0]))+90)%360,210)
        self.assertAlmostEqual(meta['predicted_moon_camera_xy'][1],50)

    def fixture(self,root):
        import tifffile
        site='WClawn';reg=root/'source';base=reg/'runs'/site;dest=base/site/'registration';dest.mkdir(parents=True)
        raw=root/'raw'/site;raw.mkdir(parents=True)
        image=np.arange(10000,dtype=np.float32).reshape(100,100)+1
        tifffile.imwrite(raw/'1.tif',image)
        (raw/'linear_raw_manifest.json').write_text('{}')
        np.savez_compressed(dest/'linear_stokes.npz',I=image,Q=.1*image,U=.15*image)
        lock=dict(site=site,status='numerically_locked_awaiting_visual_review',transforms=[{'unchanged':True}],
            reference_center_xy_px=[50,50],rim_radius_px=40,linear_raw_manifest_sha256=m.digest(raw/'linear_raw_manifest.json'),
            astronomical_anchor={'moon_elevation_deg':48.59037789,'moon_compass_azimuth_deg':210})
        path=base/'locked_site_references.json';path.write_text(json.dumps({'sites':{site:lock}}))
        (base/'reference_sha256.json').write_text(json.dumps({str(f.relative_to(base)):m.digest(f) for f in base.rglob('*') if f.is_file()}))
        summary=reg/'registration_batch_summary.json';summary.write_text(json.dumps({'sites':[dict(site=site,status=lock['status'],lock_path=str(path.relative_to(reg)),lock_sha256=m.digest(path))]}))
        review=root/'review.json';review.write_text(json.dumps({'registration_summary_sha256':m.digest(summary),'sites':[{'site':site,'alignment':'acceptable'}]}))
        overrides=root/'overrides.json';overrides.write_text(json.dumps({site:dict(moon_pixel_xy=[70,50],note='Observed',registration_lock_sha256=m.digest(path),reference_tiff_sha256=m.digest(raw/'1.tif'),linear_raw_manifest_sha256=m.digest(raw/'linear_raw_manifest.json'),image_shape_yx=[100,100])}))
        return SimpleNamespace(registration=reg,raw_root=root/'raw',visual_review=review,moon_coordinates=overrides,sites=[site],output=root/'out'),dest/'linear_stokes.npz'

    def test_unchanged_stokes_and_transforms(self):
        with tempfile.TemporaryDirectory() as d:
            a,stokes=self.fixture(Path(d));before=m.digest(stokes);m.apply(a)
            self.assertEqual(before,m.digest(a.output/'WClawn/registration/linear_stokes.npz'))
            lock=json.loads((a.output/'locked_site_references.json').read_text())['sites']['WClawn']
            self.assertEqual(lock['transforms'],[{'unchanged':True}])
            self.assertIn('camera_to_north_xy',lock)
            self.assertTrue((a.output/'WClawn/registered_dop_and_north.png').exists())

    def test_changed_raw_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as d:
            a,_=self.fixture(Path(d));(a.raw_root/'WClawn/1.tif').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'hash mismatch'):m.apply(a)
            self.assertFalse(a.output.exists())


if __name__=='__main__':unittest.main()
