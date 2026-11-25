"""Regression checks for unmapped, read-only RAW inspection."""
import contextlib
import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('raw_inspection',ROOT/'src/prepare_spider_linear_raw.py')
raw=importlib.util.module_from_spec(spec)
spec.loader.exec_module(raw)


class RawInspectionTests(unittest.TestCase):
    def test_cottage_names_sorted_by_metadata_not_museum_mapping(self):
        with tempfile.TemporaryDirectory() as td:
            folder=Path(td)
            for name in ['COTTAGE_20.NEF','COTTAGE_10.nef','COTTAGE_10.DNG']:
                (folder/name).write_bytes(b'fixture')
            def meta(p):
                return dict(datetime_original='2018:08:21 22:37:'+('47' if p.stem.endswith('20') else '58'),
                    exposure_seconds=10,iso=3200,f_number=3.5,subsecond_original='')
            with patch.object(raw,'read_metadata',side_effect=meta),patch.object(sys,'argv',
                    ['prepare','--raw-dir',td,'--inspect']),contextlib.redirect_stdout(io.StringIO()) as output:
                raw.main()
            text=output.getvalue()
            self.assertLess(text.index('COTTAGE_20.NEF'),text.index('COTTAGE_10.nef'))
            self.assertNotIn('DSC_5245',text)
            self.assertNotIn('COTTAGE_10.DNG',text)
            self.assertEqual(len(list(folder.iterdir())),3)

    def test_explicit_raw_kind_and_missing_files(self):
        with tempfile.TemporaryDirectory() as td:
            folder=Path(td)
            with self.assertRaisesRegex(ValueError,'No matching'):
                raw.inspect_available_sources(folder)
            (folder/'other.DNG').write_bytes(b'fixture')
            with patch.object(raw,'read_metadata',return_value={'datetime_original':'2018:08:21 22:37:47'}):
                self.assertEqual(raw.inspect_available_sources(folder,'dng')[0][0].name,'other.DNG')
            with self.assertRaisesRegex(ValueError,'No matching'):
                raw.inspect_available_sources(folder,'nef')


if __name__=='__main__':unittest.main()
