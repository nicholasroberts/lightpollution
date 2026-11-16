import sys,unittest,tempfile,json
from unittest.mock import patch
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from spider_natural_polarization import natural_fractional_polarization
from wang2016_sky import scattering_angle
from model_spider_nightly_prc import body_field,integrated_pattern
from plot_spider_spectral_pilot import sky_dop_snapshot

class ForwardLawTests(unittest.TestCase):
    def test_shared_phase_and_amplitude(self):
        th=np.linspace(.1,1.,31);ph=np.linspace(-2,2,31)
        old=natural_fractional_polarization(th,ph,1.1,.7,10,'berry_envelope')
        new=natural_fractional_polarization(th,ph,1.1,.7,10)
        gamma=scattering_angle(th,ph,1.1,.7)
        np.testing.assert_allclose(new,old/(1+np.cos(gamma)**2))
        np.testing.assert_allclose(np.angle(new),np.angle(old))
    def test_snapshot_reads_old_provenance_as_legacy(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'s.csv';p.write_text('wavelength_nm,natural,alan,sensitivity\n400,1,1,1\n500,1,1,1\n600,1,1,1\n')
            params=dict(brightness_model='constant',neutral_separation_deg=10)
            provenance=dict(parameters=params,sites=[dict(site='Cottage',anchor=dict(moon_elevation_deg=20,moon_compass_azimuth_deg=120))])
            legacy,geometry=sky_dop_snapshot(provenance,p,'Cottage',1,0,grid=31)
            self.assertEqual(geometry['natural_dop_law'],'berry_envelope')
            params['natural_dop_law']='rayleigh_normalized'
            normalized,geometry=sky_dop_snapshot(provenance,p,'Cottage',1,0,grid=31)
            good=np.isfinite(legacy[0])&np.isfinite(normalized[0])
            self.assertTrue(np.all(normalized[0][good]<=legacy[0][good]+1e-12))
    def test_default_field_matches_explicit_normalized(self):
        rays,axis=body_field(30,0,128)
        self.assertAlmostEqual(integrated_pattern(rays,axis,25,90,10),integrated_pattern(rays,axis,25,90,10,natural_dop_law='rayleigh_normalized'))
    def test_runner_reuses_alan_bytes_without_adjustment(self):
        import run_independent_site_alan as runner
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);(root/'batch').mkdir();configs=root/'configs';configs.mkdir()
            content='{ "historical_target": 0.45 }\n'
            for name in ['HPS_site_alan.json','LED_mid_site_alan.json']:(configs/name).write_text(content)
            args=['runner','--batch',str(root/'batch'),'--spectra',str(root/'spectra'),'--moon-anchors',str(root/'anchors.json'),'--alan-configs',str(configs),'--levels','0','5','--output',str(root/'output')]
            with patch.object(sys,'argv',args),patch.object(runner.subprocess,'run') as run:
                runner.main()
            for name in ['HPS_site_alan.json','LED_mid_site_alan.json']:self.assertEqual((root/'output/alan_configs'/name).read_text(),content)
            commands=[c.args[0] for c in run.call_args_list]
            self.assertFalse(any('adjust_cottage_alan.py' in str(c) or 'prepare_site_alan_scenarios.py' in str(c) for c in commands))
            for command in commands[:2]:self.assertIn('rayleigh_normalized',command)
            notes=json.loads((root/'output/scenario_notes.json').read_text());self.assertEqual(notes['status'],'FROZEN_PREVIOUS_ALAN_SCENARIO')
    def test_invalid_law_rejected(self):
        with self.assertRaises(ValueError):natural_fractional_polarization(.3,.4,1.,.2,10,'unknown')

if __name__=='__main__':unittest.main()
