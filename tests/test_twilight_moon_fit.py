import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from fit_spider_twilight_moon import component,mixture,fit_model,errors
from audit_spider_lunar_aop import rays,fixed_moon_in_camera

class MixtureTest(unittest.TestCase):
    def test_recovers_known_mixture_and_beats_single_sources(self):
        y,x=np.indices((60,60));view=rays(np.stack([x,y],-1),[30,30],40)
        rot=[[1,0,0],[0,1,0]]
        m=fixed_moon_in_camera(dict(moon_elevation_deg=41.2,moon_compass_azimuth_deg=187.15),rot)
        s=fixed_moon_in_camera(dict(moon_elevation_deg=-10.23,moon_compass_azimuth_deg=319.08),rot)
        moon,sun=component(view,m,.8),component(view,s,.8)
        intensity,z=mixture(moon,sun,1.7,.4);obs=abs(z)/intensity
        mask=np.hypot(x-30,y-30)<25;test=mask&((x//8+y//8)%3==0);train=mask&~test
        fit,I,Z=fit_model(obs,moon,sun,train,'mixed')
        np.testing.assert_allclose(fit['parameters'],[1.7,.4],atol=1e-5)
        self.assertLess(errors(obs,abs(Z)/I,test)['rmse'],1e-7)
        for kind in ['moon','sun']:
            _,I,Z=fit_model(obs,moon,sun,train,kind)
            self.assertGreater(errors(obs,abs(Z)/I,test)['rmse'],.01)
    def test_stokes_addition_cancels_orthogonal_components(self):
        I,z=mixture((np.array([1.]),np.array([.4+0j])),(np.array([1.]),np.array([-.4+0j])),1,0)
        np.testing.assert_allclose(I,[2]);np.testing.assert_allclose(z,[0])

if __name__=='__main__':unittest.main()
