import sys,unittest,tempfile,json
from unittest.mock import patch
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from run_spider_heading_convergence import grid_statistics,model_command
from extend_spider_nightly_sites import extend_scales,site_timezone,prepare_anchors

class HeadingExtensionTests(unittest.TestCase):
    def test_nested_grids_match_independent_selection(self):
        h=np.arange(360);v=np.sin(np.deg2rad(2*(h-3)))
        stats=grid_statistics(h,v)
        for step in [10,5,1]:self.assertAlmostEqual(stats[step]['median_abs_prc'],np.median(abs(v[::step])))
        self.assertNotEqual(stats[10]['median_abs_prc'],stats[1]['median_abs_prc'])
    def test_incomplete_grid_rejected(self):
        with self.assertRaises(ValueError):grid_statistics(np.arange(359),np.ones(359))
    def test_preserves_settings_in_fine_command(self):
        p=dict(spectra='HPS.csv',samples=1024,heading_step=10,step_minutes=20,alan_ratios=[0,5,6,7],natural_dop_law='rayleigh_normalized')
        command=model_command(p,'output',1)
        self.assertEqual(command[command.index('--heading-step')+1],'1')
        self.assertEqual(command[command.index('--samples')+1],'1024')
        self.assertEqual(p['heading_step'],10)
    def test_new_scales_use_radiance_and_preserve_existing(self):
        cfg=dict(alan_scale_by_site={'Cottage':.2,'oxfordbussinesspark':1.},spectral_band_factor=1.6)
        refs=[dict(site_id='s02',source_filename='oxfordbussinesspark.csv',original_multiyear_monthly_radiance_mean='40',original_multiyear_monthly_radiance_sd='1'),dict(site_id='s14',source_filename='Cottage.csv',original_multiyear_monthly_radiance_mean='.4',original_multiyear_monthly_radiance_sd='.1'),dict(site_id='s17',source_filename='slabhouse.csv',original_multiyear_monthly_radiance_mean='.5',original_multiyear_monthly_radiance_sd='.1')]
        out,inventory=extend_scales(cfg,refs)
        self.assertEqual(out['alan_scale_by_site']['Cottage'],.2)
        self.assertAlmostEqual(out['alan_scale_by_site']['slabhouse'],.5/40)
        self.assertEqual(out['spectral_band_factor'],1.6)
        self.assertNotIn('slabhouse',cfg['alan_scale_by_site'])
    def test_export_timestamp_stays_pending(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);folder=root/'slabhouse/registration';folder.mkdir(parents=True)
            (folder/'rotation_analysis.json').write_text(json.dumps({'inputs':[{'path':'/unused/1.tif'}]}))
            refs=[dict(source_filename='slabhouse.csv',site_id='s17',field_date='2019-05-21',latitude='51',longitude='-2')]
            with patch('anchor_moon_from_local_time.capture_time',return_value=('2019:05:21 00:49:10','TIFF DateTime fallback: may be export metadata')):
                anchors,status=prepare_anchors(root,refs,{})
            self.assertEqual(anchors,{})
            self.assertIn('fallback',status[0]['reason'])
    def test_verified_exif_gets_forward_only_ephemeris(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);folder=root/'slabhouse/registration';folder.mkdir(parents=True)
            (folder/'rotation_analysis.json').write_text(json.dumps({'inputs':[{'path':'/unused/1.tif'}]}))
            refs=[dict(source_filename='slabhouse.csv',site_id='s17',field_date='2019-05-21',latitude='51',longitude='-2')]
            with patch('anchor_moon_from_local_time.capture_time',return_value=('2019:05:21 00:49:10','EXIF DateTimeOriginal')),patch('wang2016_sky.lunar_track',return_value=(np.array([11.]),np.array([149.]))):
                anchors,status=prepare_anchors(root,refs,{})
            self.assertEqual(anchors['slabhouse']['utc'],'2019-05-20T23:49:10')
            self.assertIsNone(anchors['slabhouse']['candidate_pixel_xy'])
    def test_timezones(self):
        self.assertEqual(site_timezone('s10'),'Europe/Stockholm')
        self.assertEqual(site_timezone('s12'),'Europe/Paris')
        self.assertEqual(site_timezone('s17'),'Europe/London')

if __name__=='__main__':unittest.main()
