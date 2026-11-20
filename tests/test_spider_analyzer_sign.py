import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from correct_spider_analyzer_sign import evidence,apply,digest


def lock_for(angles):
    transforms=[]
    for n,angle in enumerate(angles,1):
        c,s=np.cos(np.deg2rad(angle)),np.sin(np.deg2rad(angle))
        transforms.append(dict(filename=str(n)+'.tif',source_to_reference_xy=[[c,-s,0],[s,c,0]],quality={'numerical_pass':True}))
    return dict(status='numerically_locked_awaiting_visual_review',analyzer_groups=[[1,2,3,4]],transforms=transforms)


class AnalyzerSignTests(unittest.TestCase):
    def test_geometry_and_linear_polarizer_recovery(self):
        for sign in [1,-1]:
            angles=sign*np.array([0.,45.,90.,135.]);lock=lock_for(angles)
            inferred,_=evidence(lock);self.assertEqual(inferred,sign)
            intensity,q,u=1.,.2,.3
            data=.5*(intensity+q*np.cos(2*np.deg2rad(angles))+u*np.sin(2*np.deg2rad(angles)))
            stored=np.array([data[0]+data[2],data[0]-data[2],data[1]-data[3]])
            np.testing.assert_allclose([stored[0],stored[1],inferred*stored[2]],[intensity,q,u])
            for angle,t in zip(angles,lock['transforms']):
                r=np.asarray(t['source_to_reference_xy'])[:,:2]
                ref=np.array([.3,.7]);moving=r.T@ref
                np.testing.assert_allclose(r@moving,ref)
                self.assertAlmostEqual(np.degrees(np.arctan2((r@[1,0])[1],(r@[1,0])[0])),angle)

    def test_real_site_rotation_sequences(self):
        for angles in [[0,-44.55568046,-88.45974322,-138.36277284],[0,-45.77839688,-89.40346095,-133.30238023]]:
            sign,trace=evidence(lock_for(angles));self.assertEqual(sign,-1)
            self.assertLess(trace[0]['maximum_nominal_angle_deviation_deg'],4)

    def test_ambiguous_and_double_correction_rejected(self):
        with self.assertRaises(ValueError):evidence(lock_for([0,10,20,30]))
        lock=lock_for([0,-45,-90,-135]);lock['analyzer_stokes_convention']={'applied_U_multiplier':-1}
        with self.assertRaises(ValueError):evidence(lock)

    def test_derived_files_preserve_source_I_Q_and_dolp(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);reg=root/'site/registration';reg.mkdir(parents=True)
            i=np.ones((4,4));q=.2*i;u=.3*i
            np.savez_compressed(reg/'linear_stokes.npz',I=i,Q=q,U=u)
            lock=lock_for([0,-45,-90,-135]);source=root/'locked_site_references.json'
            source.write_text(json.dumps({'sites':{'site':lock}}))
            (root/'reference_sha256.json').write_text(json.dumps({str(f.relative_to(root)):digest(f) for f in root.rglob('*') if f.is_file()}))
            before=digest(reg/'linear_stokes.npz')
            apply(SimpleNamespace(reference_lock=source,sites=['site'],output=root/'out'))
            self.assertEqual(digest(reg/'linear_stokes.npz'),before)
            with np.load(root/'out/site/registration/linear_stokes.npz') as z:
                np.testing.assert_array_equal(z['I'],i);np.testing.assert_array_equal(z['Q'],q);np.testing.assert_array_equal(z['U'],-u)
                np.testing.assert_array_equal(np.hypot(z['Q'],z['U'])/z['I'],np.hypot(q,u)/i)
            corrected=json.loads((root/'out/locked_site_references.json').read_text())['sites']['site']
            self.assertEqual(corrected['transforms'],lock['transforms'])


if __name__=='__main__':unittest.main()
