"""Reference polygons filter registration matches without changing image coordinates."""
import ast
from pathlib import Path
import unittest
import numpy as np

path=Path(__file__).resolve().parents[1]/'src/lock_spider_site_registration.py'
node=next(n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='filter_static_regions')
namespace={'np':np}
exec(compile(ast.Module(body=[node],type_ignores=[]),str(path),'exec'),namespace)
filter_regions=namespace['filter_static_regions']


class StaticRegionTests(unittest.TestCase):
    def test_reference_selection_preserves_correspondences(self):
        target=np.vstack([np.c_[np.arange(20)+1,np.ones(20)],[[100,100]]])
        source=target+2
        a,b=filter_regions(source,target,dict(regions_xy=[[[0,0],[30,0],[30,10],[0,10]]]))
        self.assertEqual(len(a),20)
        np.testing.assert_allclose(a-b,2)

    def test_insufficient_or_invalid_support_rejected(self):
        p=np.c_[np.arange(20)+1,np.ones(20)]
        for regions in [[],[[[0,0],[1,1]]],[[[0,0],[2,0],[2,2]]]]:
            with self.assertRaises(ValueError):filter_regions(p,p,dict(regions_xy=regions))


if __name__=='__main__':unittest.main()
