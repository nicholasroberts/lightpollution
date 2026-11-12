import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from batch_spider_view_directions import overrides,extract_rows

class DirectionBatchTests(unittest.TestCase):
    def test_override_only_known_sites(self):
        self.assertIn('Gavarnie',overrides(['Gavarnie=/tmp/example']))
        with self.assertRaises(ValueError):overrides(['wells=/tmp/example'])
    def test_extract_keeps_historical_and_band_separate(self):
        summary=dict(site_id='s12',selected_band=dict(elevation_deg=45,azimuth_deg=120),
            zenith_and_band_results=[dict(field='band_elevation_scan',method='linear_stokes',
            response_mode='stokes_intensity',prc=.603,elevation_deg=45,azimuth_deg=120,
            mean_dop=.411,physical_fraction=1)])
        row=extract_rows('Gavarnie',summary,.470,Path('/tmp/example'))[0]
        self.assertEqual(row['historical_prc'],.470)
        self.assertEqual(row['prc'],.603)
        self.assertEqual(row['field'],'band')
        self.assertAlmostEqual(row['difference_from_historical'],.133)

if __name__=='__main__':unittest.main()
