import sys,json,tempfile,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from plot_spider_site_details import measured_camera_frame,selected_rows

class SiteDetailsTest(unittest.TestCase):
    def test_camera_dop_without_invented_orientation(self):
        with tempfile.TemporaryDirectory() as root:
            folder=Path(root)/'site/registration';folder.mkdir(parents=True)
            (folder/'rotation_analysis.json').write_text(json.dumps(dict(reference_center_xy_px=[50,50],radius_60_px=30)))
            i=np.ones((101,101));np.savez(folder/'linear_stokes.npz',I=i,Q=i*.3,U=i*.4)
            dop,meta=measured_camera_frame(Path(root),'site',grid=41)
            np.testing.assert_allclose(dop[np.isfinite(dop)],.5)
            self.assertIn('north unknown',meta['orientation'])
            self.assertTrue(np.isnan(dop[0,0]))

    def test_scenario_selection_and_time_order(self):
        base=dict(site='site',atmospheric_factor='.8',alan_band_ratio='5',body_tilt_deg='0')
        rows=[dict(base,utc='2020-01-01T23:00:00'),dict(base,utc='2020-01-01T22:00:00'),dict(base,utc='2020-01-01T21:00:00',alan_band_ratio='6')]
        result=selected_rows(rows,'site',.8,5)
        self.assertEqual(len(result),2)
        self.assertEqual(result[0]['utc'],'2020-01-01T22:00:00')

if __name__=='__main__':unittest.main()
