import sys,unittest,tempfile
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from prepare_spider_linear_raw import blue_bayer_plane,output_shape,select_sources,MUSEUM_MAPPING


class RawBlueTest(unittest.TestCase):
    def test_blue_extraction_is_linear_and_preserves_negative_counts(self):
        colors=np.array([[0,1],[3,2]]);black=[64,65,66,67]
        data=np.full((6,8),1000,dtype=np.uint16);data[1::2,1::2]=np.array([[60,100,200,300]]*3)
        out,meta=blue_bayer_plane(data,colors,black,b'RGBG')
        np.testing.assert_array_equal(out[0],[-6,34,134,234]);self.assertEqual(meta['blue_offset_yx'],[1,1])
        self.assertEqual(out.dtype,np.float32)

    def test_known_polarization_survives_black_subtraction_without_wb(self):
        angles=np.deg2rad([0,45,90,135]);I,Q,U=1200.,240.,360.;black=[100,110,120,130];samples=[]
        for angle in angles:
            # Float mosaic permits an exact arithmetic regression independent of quantization.
            data=np.ones((8,8))*200;data[1::2,0::2]=black[2]+.5*(I+Q*np.cos(2*angle)+U*np.sin(2*angle))
            out,_=blue_bayer_plane(data,np.array([[1,0],[2,3]]),black,b'RGBG');samples.append(out[0,0])
        measured=np.array([samples[0]+samples[2],samples[0]-samples[2],samples[1]-samples[3]])
        np.testing.assert_allclose(measured,[I,Q,U],atol=1e-5)

    def test_exact_mapping_and_preference(self):
        with tempfile.TemporaryDirectory() as root:
            folder=Path(root)
            for suffix in ['.NEF','.dng']:(folder/('DSC_5245'+suffix)).touch()
            pair=select_sources(folder,MUSEUM_MAPPING[:1])[0]
            self.assertEqual(pair[1].suffix,'.NEF')
            with self.assertRaises(ValueError):select_sources(folder,MUSEUM_MAPPING)
        self.assertEqual(output_shape((2752,4128),2000),(1333,2000))


if __name__=='__main__':unittest.main()
