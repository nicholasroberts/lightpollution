"""Saved Oxford inputs must retain their signed Stokes and provenance."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np

spec=importlib.util.spec_from_file_location('runner',Path(__file__).resolve().parents[1]/'src/run_paper1_oxford_natural_fits.py')
runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)

class ReuseTests(unittest.TestCase):
    def fixture(self,root):
        registration=root/'registration';registration.mkdir()
        locks={}
        for site in runner.SITES:
            locks[site]=dict(radiometry_status='black_subtracted_linear_camera_green',
                analyzer_stokes_convention={'applied_U_multiplier':-1},
                moon_pixel_xy=[10,20],astronomical_anchor={'moon_compass_azimuth_deg':180})
            path=registration/site/'registration/linear_stokes.npz'
            path.parent.mkdir(parents=True)
            np.savez(path,I=np.ones((3,3)),Q=np.full((3,3),.2),U=np.full((3,3),-.3))
        lock=registration/'locked_site_references.json'
        lock.write_text(json.dumps({'sites':locks}))
        checks={str(p.relative_to(registration)):runner.digest(p) for p in registration.rglob('*') if p.is_file()}
        (registration/'reference_sha256.json').write_text(json.dumps(checks))
        for site in runner.SITES:
            p=root/site/'fit/twilight_moon_fit.json';p.parent.mkdir(parents=True)
            p.write_text(json.dumps(dict(source_lock_sha256=runner.digest(lock),
                stokes_sha256=runner.digest(registration/site/'registration/linear_stokes.npz'),
                camera_channel='green',handedness_hypothesis='reversed')))
        return lock
    def test_signed_U_is_copied_without_second_reversal(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);self.fixture(root)
            lock,prepared=runner.prepare_inputs(root)
            site,path,_,_,arrays=prepared[0]
            runner.export_reused_stokes(root/'export',site,path,lock,arrays)
            with np.load(root/'export/measured_stokes_and_polarization.npz') as z:
                for key,out in [('I','I_pair'),('Q','Q'),('U','U')]:
                    np.testing.assert_array_equal(arrays[key],z[out])
                self.assertTrue(np.all(z['U']<0))
    def test_changed_stokes_rejected_before_export(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);self.fixture(root)
            path=root/'registration'/runner.SITES[0]/'registration/linear_stokes.npz'
            path.write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'checksum'):
                runner.prepare_inputs(root)
    def test_unlinked_pilot_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);self.fixture(root)
            path=root/runner.SITES[0]/'fit/twilight_moon_fit.json'
            data=json.loads(path.read_text());data['source_lock_sha256']='other'
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError,'does not match'):
                runner.prepare_inputs(root)

if __name__=='__main__':unittest.main()
