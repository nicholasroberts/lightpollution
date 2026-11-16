import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from prepare_spider_linear_raw import green_bayer_plane,blue_bayer_plane
class GreenTest(unittest.TestCase):
    def test_individual_black_subtraction_and_green_average(self):
        colors=np.array([[0,1],[3,2]]);m=np.tile([[100,120],[150,200]],(4,5))
        green,meta=green_bayer_plane(m,colors,[10,20,30,50],b'RGBG')
        np.testing.assert_allclose(green,100)
        self.assertEqual(green.shape,(4,5));self.assertEqual(meta['green_pair_median_fractional_difference'],0)
        blue,_=blue_bayer_plane(m,colors,[10,20,30,50],b'RGBG');np.testing.assert_allclose(blue,170)
    def test_phase_and_unequal_greens(self):
        colors=np.array([[3,2],[0,1]]);m=np.tile([[70,100],[100,110]],(3,4))
        green,meta=green_bayer_plane(m,colors,[0,10,0,20],b'RGBG')
        np.testing.assert_allclose(green,75);self.assertAlmostEqual(meta['green_pair_median_fractional_difference'],-2/3)
if __name__=='__main__':unittest.main()
