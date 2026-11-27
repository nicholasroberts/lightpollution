"""Figure 2 must use the full approved 400–700 nm energy band."""
import ast
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np

source=Path(__file__).resolve().parents[1]/'src/plot_paper1_fig02a.py'
scope=dict(np=np,json=json,digest=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest())
tree=ast.parse(source.read_text())
nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ('integral','load_spectrum','spectral_inputs')]
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),scope)

class Figure2SpectralBandTests(unittest.TestCase):
    def inputs(self,root,lo):
        wave=np.arange(float(lo),751.)
        for name,peak in [('HPS',590.),('LED_mid',450.)]:
            np.savetxt(root/(name+'_visible_input.csv'),
                       np.column_stack([wave,np.ones_like(wave),np.exp(-((wave-peak)/30.)**2)]),
                       delimiter=',',header='wavelength_nm,natural,alan',comments='')
        (root/'spectral_provenance.json').write_text('{}')

    def test_full_band_and_unit_energy(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.inputs(root,352)
            wave,spectra,report=scope['spectral_inputs'](root)
            self.assertEqual([wave[0],wave[-1]],[400.,700.])
            self.assertEqual(report['spectral_band_nm'],[400.,700.])
            for spectrum in spectra.values():
                self.assertAlmostEqual(scope['integral'](spectrum,wave),1.,places=12)

    def test_old_450_inputs_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.inputs(root,450)
            with self.assertRaisesRegex(ValueError,'full 400–700'):
                scope['spectral_inputs'](root)

if __name__=='__main__':unittest.main()
