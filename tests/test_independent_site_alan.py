import sys,unittest,tempfile,json
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from prepare_site_alan_scenarios import overlap_conversion
from plot_spider_spectral_pilot import measured_north_up
from reconstruct_spider_eye_contrast import digest

class IndependentAlanTests(unittest.TestCase):
    def test_source_scale_cancels_dnb_conversion(self):
        source=np.column_stack([np.arange(300.,901.),np.ones(601)])
        response=np.array([[300.,1.],[900.,1.]])
        factor,_=overlap_conversion(source,response,0)
        scaled=source.copy();scaled[:,1]*=7
        self.assertAlmostEqual(factor,.5)
        self.assertAlmostEqual(overlap_conversion(scaled,response,0)[0],factor)
    def test_spectral_redistribution_changes_conversion(self):
        w=np.arange(300.,901.);rsr=np.array([[300.,.1],[900.,1.]])
        blue=np.column_stack([w,np.exp(-((w-460)/25)**2)])
        red=np.column_stack([w,np.exp(-((w-650)/25)**2)])
        self.assertGreater(overlap_conversion(blue,rsr,2)[0],overlap_conversion(red,rsr,2)[0])
    def test_measured_reprojection_preserves_uniform_dop(self):
        with tempfile.TemporaryDirectory() as td:
            folder=Path(td)/'Cottage/registration';folder.mkdir(parents=True)
            report=folder/'rotation_analysis.json';report.write_text(json.dumps({'reference_center_xy_px':[25,25],'radius_60_px':20/np.sqrt(2)}))
            i=np.ones((51,51));np.savez(folder/'linear_stokes.npz',I=i,Q=.2*i,U=0*i)
            geometry={'anchor':{'source_report_sha256':digest(report),'moon_compass_azimuth_deg':90,'candidate_pixel_xy':[40,25]}}
            image,_=measured_north_up(Path(td),'Cottage',geometry,grid=41)
            np.testing.assert_allclose(image[np.isfinite(image)],.2)
            geometry['anchor']['source_report_sha256']='different'
            with self.assertRaises(ValueError):measured_north_up(Path(td),'Cottage',geometry,grid=41)

if __name__=='__main__':unittest.main()
