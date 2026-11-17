import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from prepare_spider_raw_batch import make_plans

def rows(folder,names,seconds='99/10'):
    return [dict(path=f'data/raw/skyimages/{folder}/{name}.NEF',size_bytes=10,metadata={
        'EXIF DateTimeOriginal':f'2019:05:14 22:{i:02d}:00','EXIF OffsetTimeOriginal':'+01:00',
        'Image Model':'NIKON D850','EXIF ExposureTime':seconds,'EXIF ISOSpeedRatings':'3200','EXIF FNumber':'7/2'}) for i,name in enumerate(names)]
class BatchTest(unittest.TestCase):
    def test_retained_order_ignores_filename_gaps(self):
        r=rows('WC_pitch',['DSC_'+str(n) for n in [1,2,3,7,8,9,12,13,14,20,21,22]])
        p=make_plans(r,Path('/tmp/raw'))[0]
        self.assertEqual(p['site'],'WCpitch');self.assertEqual(p['mapping'][3]['stem'],'DSC_7')
        self.assertEqual(p['analyzer_groups'],[[1,4,7,10],[2,5,8,11],[3,6,9,12]])
    def test_waitrose_single_quartet(self):
        p=make_plans(rows('Waitrose',['a','b','c','d']),Path('/tmp/raw'))[0]
        self.assertEqual(p['analyzer_groups'],[[1,2,3,4]])
    def test_known_anomaly_excluded_without_renumbering(self):
        r=rows('slabhouse',['DSC_'+str(n) for n in range(5395,5407)],'249/10');r[9]['metadata']['EXIF ExposureTime']='49/2'
        p=make_plans(r,Path('/tmp/raw'))[0]
        self.assertEqual(p['excluded_exposures'],[10]);self.assertEqual(p['mapping'][10]['exposure'],11)
        self.assertEqual(len(p['analyzer_groups']),2)
    def test_unknown_changes_and_reversed_time_rejected(self):
        r=rows('test',['a','b','c','d']);r[2]['metadata']['EXIF ExposureTime']='1'
        with self.assertRaises(ValueError):make_plans(r,Path('/tmp/raw'))
        r=rows('test',['a','b','c','d']);r[2]['metadata']['EXIF DateTimeOriginal']=r[0]['metadata']['EXIF DateTimeOriginal']
        with self.assertRaises(ValueError):make_plans(r,Path('/tmp/raw'))
if __name__=='__main__':unittest.main()
