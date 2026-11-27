"""Checks for independent spatial lamp fits; synthetic data only."""
import ast
from pathlib import Path
import unittest
import numpy as np
from scipy.optimize import brentq,least_squares

# Load the numerical routine without importing plotting/archive dependencies.
source=Path(__file__).resolve().parents[1]/'src/fit_paper1_fig02b.py'
node=next(n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='fit_amount')
scope=dict(np=np,least_squares=least_squares)
exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),scope)
fit_amount=scope['fit_amount']

class SpatialLampFitTests(unittest.TestCase):
    def test_known_amount_and_independent_spectral_alternatives(self):
        i=np.linspace(1.,2.,100);p=np.linspace(.3,.9,100)
        observed=p/(i+.7*2.)
        hps=fit_amount(i,p,observed,2.)
        led=fit_amount(i,p,observed,3.)
        self.assertAlmostEqual(hps['amount'],.7,places=7)
        self.assertAlmostEqual(led['amount'],.7*2/3,places=7)
        self.assertLess(hps['objective_mean_squared_DoLP_residual'],1e-15)

    def test_spatial_objective_improves_on_mean_matching(self):
        i=np.ones(100);p=np.linspace(.1,.9,100)
        observed=np.full(100,.25)
        mean_amount=brentq(lambda a:np.mean(p/(i+a))-.25,0,10)
        fitted=fit_amount(i,p,observed,1.)
        old_cost=np.mean((p/(i+mean_amount)-observed)**2)
        self.assertLess(fitted['objective_mean_squared_DoLP_residual'],old_cost)
        self.assertGreater(abs(fitted['achieved_mean_dolp']-.25),.01)

    def test_zero_boundary_when_addition_cannot_help(self):
        fitted=fit_amount(np.ones(40),np.full(40,.3),np.full(40,.4),2.)
        self.assertEqual(fitted['amount'],0.)
        self.assertTrue(fitted['boundary_zero'])

if __name__=='__main__':unittest.main()
