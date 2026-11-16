import sys,unittest,tempfile,json
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from audit_spider_lunar_aop import change_residual,lunar_phase,fixed_moon_in_camera,block_mean,main,digest
from spider_fisheye_registration import rays,pixels


class AoPGeometryTest(unittest.TestCase):
    def test_global_analyzer_zero_cancels_across_wrapped_angles(self):
        y,x=np.indices((30,40));angle=.08*x+.03*y
        predicted=np.exp(2j*angle);measured=np.exp(2j*(angle+1.3));mask=np.ones(angle.shape,bool)
        for axis in [0,1]:
            residual=change_residual(measured,predicted,mask,axis,4)
            np.testing.assert_allclose(residual,0,atol=1e-14)
        reversed_residual=change_residual(np.conj(measured),predicted,mask,1,4)
        self.assertGreater(np.median(abs(reversed_residual)),.1)

    def test_ray_basis_at_zenith_and_fixed_moon_transform(self):
        moon=fixed_moon_in_camera({'moon_elevation_deg':45,'moon_compass_azimuth_deg':180},[[1,0,0],[0,1,0]])
        self.assertGreater(moon[1],0);self.assertAlmostEqual(np.linalg.norm(moon),1)
        np.testing.assert_allclose(lunar_phase(np.array([[0.,0.,1.]]),moon),[1+0j],atol=1e-12)
        np.testing.assert_allclose(block_mean(np.arange(16.).reshape(4,4),2),[[2.5,4.5],[10.5,12.5]])

    def test_checksum_locked_synthetic_cli_preserves_unknown_zero(self):
        with tempfile.TemporaryDirectory() as root:
            root=Path(root);reg=root/'museum/registration';reg.mkdir(parents=True)
            c=np.array([200.,200.]);r=180.;y,x=np.indices((400,400));v=rays(np.stack([x,y],axis=-1),c,r)
            anchor={'moon_elevation_deg':45.,'moon_compass_azimuth_deg':180.,'site_id':'s15','utc':'2019-05-14T21:09:51+00:00'};moon=fixed_moon_in_camera(anchor,[[1,0,0],[0,1,0]])
            phase=lunar_phase(v,moon)*np.exp(2j*.7);I=np.ones((400,400))*100
            np.savez_compressed(reg/'linear_stokes.npz',I=I,Q=40*phase.real,U=40*phase.imag)
            site=dict(status='numerically_locked_awaiting_visual_review',radiometry_status='black_subtracted_linear_camera_blue',
                reference_center_xy_px=c.tolist(),rim_radius_px=r,eye_radius_px=r*np.sin(np.pi/12)/np.sin(np.pi/4),
                moon_pixel_xy=pixels(moon[None,:],c,r)[0].tolist(),camera_to_north_xy=[[1,0,0],[0,1,0]],astronomical_anchor=anchor)
            lock=root/'locked_site_references.json';lock.write_text(json.dumps({'sites':{'museum':site}}))
            (root/'reference_sha256.json').write_text(json.dumps({str(p.relative_to(root)):digest(p) for p in [lock,reg/'linear_stokes.npz']}))
            reference=root/'reference.csv';reference.write_text('source_filename,site_id,latitude,longitude\nmuseum.csv,s15,51.758527,-1.256267\n')
            with patch.object(sys,'argv',['audit','--reference-lock',str(lock),'--output',str(root/'audit'),'--bin-size','4','--reference',str(reference)]):main()
            result=json.loads((root/'audit/aop_geometry_summary.json').read_text())
            self.assertIsNone(result['fitted_analyzer_zero_deg']);self.assertIsNone(result['fitted_pattern_rotation_deg'])
            stat=result['metrics']['stored_U']['fixed_ephemeris_zenith_axis']['16']
            self.assertGreater(stat['pair_count'],100);self.assertLess(stat['median_abs_change_difference_deg'],.1)
            self.assertTrue((root/'audit/aop_geometry_audit.png').is_file())
            self.assertTrue((root/'audit/aop_moon_sun_control.png').is_file())
            self.assertAlmostEqual(result['solar_control']['sun_elevation_deg'],-10.2274,places=3)


if __name__=='__main__':unittest.main()
