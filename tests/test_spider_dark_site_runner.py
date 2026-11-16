import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from run_spider_dark_site_validation import corrected_anchors,rotate_point,profile_fwhm
import numpy as np

class BearingTests(unittest.TestCase):
    def test_only_bearing_changes(self):
        source={'sites':{'Skillinge':dict(candidate_pixel_xy=[637,193],moon_elevation_deg=12.5925,utc='2018-08-21T20:09:22',source_report_sha256='abc',flags=['other flag'])}}
        override=dict(source='visual overlay',status='APPROXIMATE',sites={'Skillinge':dict(candidate_pixel_xy=[650,1135],approximate_pixel_uncertainty=25)})
        out=corrected_anchors([source],override)['sites']['Skillinge']
        self.assertEqual(out['candidate_pixel_xy'],[650,1135])
        self.assertEqual(out['moon_elevation_deg'],12.5925)
        self.assertEqual(out['source_report_sha256'],'abc')
        self.assertEqual(out['utc'],source['sites']['Skillinge']['utc'])
        self.assertEqual(source['sites']['Skillinge']['candidate_pixel_xy'],[637,193])
        self.assertIn('other flag',out['flags'])
    def test_azimuth_reversal_uses_camera_centre(self):
        np.testing.assert_allclose(rotate_point([650,1135],[1000,650],180),[1350,165],atol=1e-10)
    def test_fwhm_scaling_invariant(self):
        x=np.arange(0,50,10);y=np.array([0,.5,1,.5,0])
        for scale in [1,.4]:self.assertAlmostEqual(profile_fwhm(x,y*scale)['fwhm_deg'],20)
    def test_fwhm_censoring_not_extrapolated(self):
        result=profile_fwhm([10,20,30,40],[.1,.5,1,.8])
        self.assertIsNone(result['fwhm_deg'])
        self.assertIsNone(result['right_crossing_deg'])
    def test_missing_ephemeris_rejected(self):
        with self.assertRaises(ValueError):corrected_anchors([{'sites':{}}],dict(source='x',status='x',sites={'s':dict(candidate_pixel_xy=[1,2])}))
    def test_conflicting_anchor_rejected(self):
        with self.assertRaises(ValueError):corrected_anchors([{'sites':{'s':{'x':1}}},{'sites':{'s':{'x':2}}}],dict(sites={}))

if __name__=='__main__':unittest.main()
