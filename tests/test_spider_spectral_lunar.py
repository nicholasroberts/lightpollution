import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from spider_spectral_lunar import pigment_proxy,ks91_brightness

class SpectralLunarTests(unittest.TestCase):
    def test_pigment_peak_and_uv_scenario(self):
        wave=np.arange(300.,751.)
        no= pigment_proxy(wave,uv_amplitude=0.)
        self.assertEqual(wave[np.argmax(no)],510)
        self.assertAlmostEqual(float(pigment_proxy(510)),1)
        self.assertGreater(float(pigment_proxy(365)),float(pigment_proxy(365,uv_amplitude=0)))
    def test_ks91_published_implementation_example(self):
        b=ks91_brightness(np.deg2rad(20),20,np.deg2rad(50),45,.15)
        magnitude=(20.7233-np.log(b/34.08))/.92104
        self.assertAlmostEqual(float(magnitude),19.855,places=3)
    def test_phase_and_distance(self):
        args=(0.,30.,np.pi/2)
        full=ks91_brightness(*args,0)
        self.assertGreater(full,ks91_brightness(*args,90))
        self.assertAlmostEqual(float(ks91_brightness(*args,0,distance_km=768800)/full),.25)
    def test_circumlunar_exclusion(self):
        with self.assertRaises(ValueError):ks91_brightness(0,30,.01,45)
    def test_known_capture_scaling_and_constant_alan(self):
        from model_spider_nightly_prc import response_contrast,spectral_coefficients
        co,_=spectral_coefficients()
        self.assertAlmostEqual(response_contrast(.3,1,0,co),response_contrast(.6,1,0,co,2))
        self.assertGreater(abs(response_contrast(.6,1,1,co,2)),abs(response_contrast(.3,1,1,co)))

class SkySnapshotTests(unittest.TestCase):
    def test_mixed_stokes_and_zero_alan(self):
        import tempfile
        from plot_spider_spectral_pilot import sky_dop_snapshot
        p={'sites':[{'site':'Cottage','anchor':{'moon_elevation_deg':20.,'moon_compass_azimuth_deg':150.,'utc':'2018-08-21T20:37:47'}}],
           'parameters':{'brightness_model':'constant','neutral_separation_deg':24.}}
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'s.csv'
            path.write_text('wavelength_nm,natural,alan,sensitivity\n450,1,1,1\n750,1,1,1\n')
            arrays,_=sky_dop_snapshot(p,path,'Cottage',.8,1,grid=51)
            mask=np.isfinite(arrays[0]);np.testing.assert_allclose(arrays[1][mask],arrays[0][mask]/2)
            zero,_=sky_dop_snapshot(p,path,'Cottage',.8,0,grid=51)
            np.testing.assert_allclose(zero[0],zero[1],equal_nan=True)
            self.assertLessEqual(float(np.nanmax(arrays[0])),.8)

if __name__=='__main__':unittest.main()
