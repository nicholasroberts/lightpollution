import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from run_spider_dark_site_validation import corrected_anchors

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
    def test_missing_ephemeris_rejected(self):
        with self.assertRaises(ValueError):corrected_anchors([{'sites':{}}],dict(source='x',status='x',sites={'s':dict(candidate_pixel_xy=[1,2])}))
    def test_conflicting_anchor_rejected(self):
        with self.assertRaises(ValueError):corrected_anchors([{'sites':{'s':{'x':1}}},{'sites':{'s':{'x':2}}}],dict(sites={}))

if __name__=='__main__':unittest.main()
