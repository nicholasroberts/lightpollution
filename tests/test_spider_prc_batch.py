import sys
from pathlib import Path
import tempfile
import unittest
import numpy as np
import tifffile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from batch_spider_prc_sites import locate,transferred_geometry,SITES

class BatchTests(unittest.TestCase):
    def test_locked_site_roster(self):
        self.assertEqual(len(SITES),18)
        self.assertEqual(len({v[0] for v in SITES.values()}),18)
        self.assertFalse({'wells','sandhammeren','bristol'} & {k.lower() for k in SITES})

    def test_waitrose_single_replicate(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);folder=root/'Waitrose'/'extract';folder.mkdir(parents=True)
            for n in range(1,5):(folder/f'{n}.tif').touch()
            found,files=locate(root,['Waitrose'])
            self.assertEqual(found,folder);self.assertEqual(len(files),4)
            (folder/'01.tif').touch()
            with self.assertRaises(ValueError):locate(root,['Waitrose'])

    def test_scaled_geometry_and_crop_flag(self):
        template=dict(channel=2,records=[dict(file='1.tif',shape=[1333,2000,3],fits=[dict(
            threshold_fraction=.05,center_x_px=1000,center_y_px=666.5,rim_radius_px=650)])])
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'1.tif';tifffile.imwrite(path,np.zeros((100,150,3),dtype='uint8'))
            g,flags=transferred_geometry([path],template)
            self.assertFalse(flags)
            self.assertAlmostEqual(g['records'][0]['fits'][0]['rim_radius_px'],650*100/1333)
            tifffile.imwrite(path,np.zeros((150,100,3),dtype='uint8'))
            _,flags=transferred_geometry([path],template)
            self.assertTrue(flags)

if __name__=='__main__':unittest.main()
